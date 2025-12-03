#include "dna_sequence.h"


#define SPGIST_MAX_PREFIX_LENGTH    Max((int) (BLCKSZ - 258 * 16 - 100), 32)

#define KMER_EQUAL_STRATEGY         1   /* kmer = kmer */
#define KMER_PREFIX_STRATEGY        2   /* kmer ^@ kmer (prefix) */
#define QKMER_CONTAINS_STRATEGY     3   /* qkmer @> kmer (pattern match) */

static inline uint8_t *
KMER_DATA(Kmer *k)
{
    return k->code;
}

static inline int
KMER_LEN(Kmer *k)
{
    return k->length;
}

static inline uint8_t *
QKMER_DATA(Qkmer *q)
{
    return q->code;
}

static inline int
QKMER_LEN(Qkmer *q)
{
    return q->length;
}

static bool
qkmer_base_matches(uint8_t qbase, uint8_t kbase)
{
    return (qbase & kbase) != 0;
}

static bool
qkmer_matches_kmer(const uint8_t *qdata, int qlen, 
                   const uint8_t *kdata, int klen)
{
    int i;
    
    if (qlen != klen)
        return false;
    
    for (i = 0; i < qlen; i++)
    {
        if (!qkmer_base_matches(qdata[i], kdata[i]))
            return false;
    }
    
    return true;
}


PG_FUNCTION_INFO_V1(qkmer_contains);
Datum qkmer_contains(PG_FUNCTION_ARGS)
{
    Qkmer *qk = PG_GETARG_QKMER_P(0);
    Kmer *k = PG_GETARG_KMER_P(1);
    bool result;

	result = qkmer_matches_kmer(QKMER_DATA(qk), QKMER_LEN(qk),
                                KMER_DATA(k), KMER_LEN(k));

    PG_FREE_IF_COPY(qk, 0);
    PG_FREE_IF_COPY(k, 1);
    PG_RETURN_BOOL(result);
}


static inline int
pg_cmp_s16(int16 a, int16 b)
{
    return (int32) a - (int32) b;
}


// SP-GiST Config
PG_FUNCTION_INFO_V1(spg_kmer_config);
Datum spg_kmer_config(PG_FUNCTION_ARGS)
{
    spgConfigIn *cfgin = (spgConfigIn *) PG_GETARG_POINTER(0);
    spgConfigOut *cfg = (spgConfigOut *) PG_GETARG_POINTER(1);

    cfg->prefixType = cfgin->attType;
    cfg->labelType = INT2OID;
    cfg->leafType = cfgin->attType;
    cfg->canReturnData = true;
    cfg->longValuesOK = false;
    
    PG_RETURN_VOID();
}

static Datum
formKmerDatum(const uint8_t *data, int datalen)
{
    /* Calculate total size: header + length byte + data */
    int32       totalSize = VARHDRSZ + sizeof(uint8_t) + datalen;
    Kmer       *res = (Kmer *) palloc0(totalSize);

    SET_VARSIZE(res, totalSize);
    res->length = (uint8_t) datalen;

    if (datalen > 0 && data != NULL)
        memcpy(res->code, data, datalen);

    return PointerGetDatum(res);
}

static int
commonPrefix(const uint8_t *a, const uint8_t *b, int len_a, int len_b)
{
    int i;
    int min_len = (len_a < len_b) ? len_a : len_b;

    for (i = 0; i < min_len; i++)
    {
        if (a[i] != b[i])
            break;
    }
    return i;
}


static bool
searchChar(const Datum *nodeLabels, int nNodes, int16 c, int *i)
{
    int         StopLow = 0,
                StopHigh = nNodes;
 
    while (StopLow < StopHigh)
    {
        int         StopMiddle = (StopLow + StopHigh) >> 1;
        int16       middle = DatumGetInt16(nodeLabels[StopMiddle]);
 
        if (c < middle)
            StopHigh = StopMiddle;
        else if (c > middle)
            StopLow = StopMiddle + 1;
        else
        {
            *i = StopMiddle;
            return true;
        }
    }
 
    *i = StopHigh;
    return false;
}


// SP-GiST Choose
PG_FUNCTION_INFO_V1(spg_kmer_choose);
Datum
spg_kmer_choose(PG_FUNCTION_ARGS)
{
    spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);
    Kmer       *inKmer = DatumGetKmerP(in->datum);
    uint8_t    *inData = KMER_DATA(inKmer);
    int         inLen = KMER_LEN(inKmer);
    uint8_t    *prefixData = NULL;
    int         prefixLen = 0;
    int         commonLen = 0;
    int16       nodeChar = 0;
    int         i = 0;

    /* Check for prefix match, set nodeChar to first byte after prefix */
    if (in->hasPrefix)
    {
        Kmer *prefixKmer = DatumGetKmerP(in->prefixDatum);

        prefixData = KMER_DATA(prefixKmer);
        prefixLen = KMER_LEN(prefixKmer);

        /* * commonPrefix logic: compares byte arrays. 
         * Assumes commonPrefix(char*, char*, int, int) exists or you allow casting.
         * If using the standard pg commonPrefix, cast uint8_t* to char*.
         */
        commonLen = commonPrefix(inData + in->level,
                                 prefixData,
                                 inLen - in->level,
                                 prefixLen);

        if (commonLen == prefixLen)
        {
            if (inLen - in->level > commonLen)
                nodeChar = inData[in->level + commonLen];
            else
                nodeChar = -1;
        }
        else
        {
            /* Must split tuple because incoming value doesn't match prefix */
            out->resultType = spgSplitTuple;

            if (commonLen == 0)
            {
                out->result.splitTuple.prefixHasPrefix = false;
            }
            else
            {
                out->result.splitTuple.prefixHasPrefix = true;
                out->result.splitTuple.prefixPrefixDatum =
                    formKmerDatum(prefixData, commonLen);
            }
            out->result.splitTuple.prefixNNodes = 1;
            out->result.splitTuple.prefixNodeLabels =
                (Datum *) palloc(sizeof(Datum));
            out->result.splitTuple.prefixNodeLabels[0] =
                Int16GetDatum(prefixData[commonLen]);

            out->result.splitTuple.childNodeN = 0;

            if (prefixLen - commonLen == 1)
            {
                out->result.splitTuple.postfixHasPrefix = false;
            }
            else
            {
                out->result.splitTuple.postfixHasPrefix = true;
                out->result.splitTuple.postfixPrefixDatum =
                    formKmerDatum(prefixData + commonLen + 1,
                                  prefixLen - commonLen - 1);
            }

            PG_RETURN_VOID();
        }
    }
    else if (inLen > in->level)
    {
        nodeChar = inData[in->level];
    }
    else
    {
        nodeChar = -1;
    }

    /* Look up nodeChar in the node label array */
    if (searchChar(in->nodeLabels, in->nNodes, nodeChar, &i))
    {
        /* Descend to existing node */
        int         levelAdd;

        out->resultType = spgMatchNode;
        out->result.matchNode.nodeN = i;
        levelAdd = commonLen;
        if (nodeChar >= 0)
            levelAdd++;
        out->result.matchNode.levelAdd = levelAdd;
        
        if (inLen - in->level - levelAdd > 0)
            out->result.matchNode.restDatum =
                formKmerDatum(inData + in->level + levelAdd,
                              inLen - in->level - levelAdd);
        else
            out->result.matchNode.restDatum =
                formKmerDatum(NULL, 0);
    }
    else if (in->allTheSame)
    {
        /* Split tuple logic for allTheSame case */
        out->resultType = spgSplitTuple;
        out->result.splitTuple.prefixHasPrefix = in->hasPrefix;
        out->result.splitTuple.prefixPrefixDatum = in->prefixDatum;
        out->result.splitTuple.prefixNNodes = 1;
        out->result.splitTuple.prefixNodeLabels = (Datum *) palloc(sizeof(Datum));
        out->result.splitTuple.prefixNodeLabels[0] = Int16GetDatum(-2);
        out->result.splitTuple.childNodeN = 0;
        out->result.splitTuple.postfixHasPrefix = false;
    }
    else
    {
        /* Add a node for the not-previously-seen nodeChar value */
        out->resultType = spgAddNode;
        out->result.addNode.nodeLabel = Int16GetDatum(nodeChar);
        out->result.addNode.nodeN = i;
    }

    PG_RETURN_VOID();
}

/* qsort comparator to sort spgNodePtr structs by "c" */
static int
cmpNodePtr(const void *a, const void *b)
{
    const spgNodePtr *aa = (const spgNodePtr *) a;
    const spgNodePtr *bb = (const spgNodePtr *) b;
 
    return pg_cmp_s16(aa->c, bb->c);
}

// SP-GiST Picksplit
PG_FUNCTION_INFO_V1(spg_kmer_picksplit);
Datum
spg_kmer_picksplit(PG_FUNCTION_ARGS)
{
    spgPickSplitIn  *in  = (spgPickSplitIn *) PG_GETARG_POINTER(0);
    spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);

    Kmer *k0       = DatumGetKmerP(in->datums[0]);
    int         i,
                commonLen;
    spgNodePtr *nodes;

	commonLen = KMER_LEN(k0);
    for (i = 1; i < in->nTuples && commonLen > 0; i++)
    {
        Kmer *ki   = DatumGetKmerP(in->datums[i]);
        int       tmp = commonPrefix(KMER_DATA(k0),
                                 KMER_DATA(ki),
                                 KMER_LEN(k0),
                                 KMER_LEN(ki));
        if (tmp < commonLen)
            commonLen = tmp;
    }

    /*
     * Limit the prefix length, if necessary, to ensure that the resulting
     * inner tuple will fit on a page.
     */
    commonLen = Min(commonLen, SPGIST_MAX_PREFIX_LENGTH);

    /* Set node prefix to be that string, if it's not empty */
    if (commonLen == 0)
    {
        out->hasPrefix = false;
    }
    else
    {
        out->hasPrefix = true;
        out->prefixDatum = formKmerDatum(KMER_DATA(k0), commonLen);
    }


    /* Extract the node label (first non-common byte) from each value */
    nodes = (spgNodePtr *) palloc(sizeof(spgNodePtr) * in->nTuples);
 
    for (i = 0; i < in->nTuples; i++)
    {
        Kmer       *ki = DatumGetKmerP(in->datums[i]);
 
        if (commonLen < KMER_LEN(ki))
            nodes[i].c = ki->code[commonLen];
        else
            nodes[i].c = -1;    /* use -1 if string is all common */
        nodes[i].i = i;
        nodes[i].d = in->datums[i];
    }


    /*
     * Sort by label values so that we can group the values into nodes.  This
     * also ensures that the nodes are ordered by label value, allowing the
     * use of binary search in searchChar.
     */
    qsort(nodes, in->nTuples, sizeof(*nodes), cmpNodePtr);
 
    /* And emit results */
    out->nNodes = 0;
    out->nodeLabels = (Datum *) palloc(sizeof(Datum) * in->nTuples);
    out->mapTuplesToNodes = (int *) palloc(sizeof(int) * in->nTuples);
    out->leafTupleDatums = (Datum *) palloc(sizeof(Datum) * in->nTuples);


    for (i = 0; i < in->nTuples; i++)
    {
        Kmer       *ki = DatumGetKmerP(nodes[i].d);
        Datum       leafD;
 
        if (i == 0 || nodes[i].c != nodes[i - 1].c)
        {
            out->nodeLabels[out->nNodes] = Int16GetDatum(nodes[i].c);
            out->nNodes++;
        }
 
        if (commonLen < KMER_LEN(ki))
			// if (ki->length - commonLen - 1 > 0)
				leafD = formKmerDatum(KMER_DATA(ki) + commonLen + 1,
						  KMER_LEN(ki) - commonLen - 1);
			// else
			// 	leadD = formKmerDatum(NULL, 0)
        else
            leafD = formKmerDatum(NULL, 0);
 
        out->leafTupleDatums[nodes[i].i] = leafD;
        out->mapTuplesToNodes[nodes[i].i] = out->nNodes - 1;
    }

    PG_RETURN_VOID();
}

// SP-GiST Inner Consistent
PG_FUNCTION_INFO_V1(spg_kmer_inner_consistent);
Datum
spg_kmer_inner_consistent(PG_FUNCTION_ARGS)
{
    spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
    spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);
    Kmer       *reconstructedValue;
    Kmer       *reconstrKmer;
    int         maxReconstrLen;
    Kmer       *prefixKmer = NULL;
    int         prefixSize = 0;
    int         i;

    /* Reconstruct the kmer at this level */
    reconstructedValue = (Kmer *) DatumGetPointer(in->reconstructedValue);
    Assert(reconstructedValue == NULL ? in->level == 0 :
           KMER_LEN(reconstructedValue) == in->level);

    maxReconstrLen = in->level + 1;
    if (in->hasPrefix)
    {
        prefixKmer = DatumGetKmerP(in->prefixDatum);
        prefixSize = KMER_LEN(prefixKmer);
        maxReconstrLen += prefixSize;
    }

    reconstrKmer = palloc(VARHDRSZ + sizeof(uint8_t) + maxReconstrLen);
    SET_VARSIZE(reconstrKmer, VARHDRSZ + sizeof(uint8_t) + maxReconstrLen);
    reconstrKmer->length = maxReconstrLen;

    if (in->level)
        memcpy(KMER_DATA(reconstrKmer),
               KMER_DATA(reconstructedValue),
               in->level);
    if (prefixSize)
        memcpy(KMER_DATA(reconstrKmer) + in->level,
               KMER_DATA(prefixKmer),
               prefixSize);

    /* Scan child nodes */
    out->nodeNumbers = (int *) palloc(sizeof(int) * in->nNodes);
    out->levelAdds = (int *) palloc(sizeof(int) * in->nNodes);
    out->reconstructedValues = (Datum *) palloc(sizeof(Datum) * in->nNodes);
    out->nNodes = 0;

    for (i = 0; i < in->nNodes; i++)
    {
        int16       nodeChar = DatumGetInt16(in->nodeLabels[i]);
        int         thisLen;
        bool        res = true;
        int         j;

        /* Handle dummy node labels */
        if (nodeChar <= 0)
            thisLen = maxReconstrLen - 1;
        else
        {
            KMER_DATA(reconstrKmer)[maxReconstrLen - 1] = nodeChar;
            thisLen = maxReconstrLen;
        }

        for (j = 0; j < in->nkeys; j++)
        {
            StrategyNumber strategy = in->scankeys[j].sk_strategy;

            switch (strategy)
            {
                case KMER_EQUAL_STRATEGY:
                    {
						Kmer       *queryKmer = DatumGetKmerP(in->scankeys[j].sk_argument);
						int         queryLen = KMER_LEN(queryKmer);
						int         r;
                        r = memcmp(KMER_DATA(reconstrKmer), 
                                   KMER_DATA(queryKmer),
                                   Min(queryLen, thisLen));

                        if (r != 0 || queryLen < thisLen)
                            res = false;
                    }
                    break;

                case KMER_PREFIX_STRATEGY:
                    {
						Kmer       *queryKmer = DatumGetKmerP(in->scankeys[j].sk_argument);
						int         queryLen = KMER_LEN(queryKmer);
						int         r;
                        /* Check if reconstructed value matches prefix so far */
                        r = memcmp(KMER_DATA(reconstrKmer),
                                   KMER_DATA(queryKmer),
                                   Min(queryLen, thisLen));

                        if (r != 0)
                            res = false;
                    }
                    break;

                case QKMER_CONTAINS_STRATEGY:
                    {
						Qkmer      *queryQkmer = DatumGetQkmerP(in->scankeys[j].sk_argument);
                        int         queryLen = QKMER_LEN(queryQkmer);
                        int         k;
                        /* Check if pattern matches so far */
                        if (thisLen > queryLen)
                        {
                            res = false;
                            break;
                        }

                        for (k = 0; k < thisLen; k++)
                        {
                            if (!qkmer_base_matches(QKMER_DATA(queryQkmer)[k],
                                                    KMER_DATA(reconstrKmer)[k]))
                            {
                                res = false;
                                break;
                            }
                        }
                    }
                    break;

                default:
                    elog(ERROR, "unrecognized strategy number: %d", strategy);
                    break;
            }

            if (!res)
                break;
        }

        if (res)
        {
            out->nodeNumbers[out->nNodes] = i;
            out->levelAdds[out->nNodes] = thisLen - in->level;
            reconstrKmer->length = thisLen;
            SET_VARSIZE(reconstrKmer, VARHDRSZ + sizeof(uint8_t) + thisLen);
            out->reconstructedValues[out->nNodes] =
                datumCopy(PointerGetDatum(reconstrKmer), false, -1);
            out->nNodes++;
        }
    }

    PG_RETURN_VOID();
}

// SP-GiST Leaf Consistent
PG_FUNCTION_INFO_V1(spg_kmer_leaf_consistent);
Datum
spg_kmer_leaf_consistent(PG_FUNCTION_ARGS)
{
    spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
    spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);
    int         level = in->level;
    Kmer       *leafValue,
               *reconstrValue = NULL;
    uint8_t    *fullValue;
    int         fullLen;
    bool        res;
    int         j;

    out->recheck = false;

    leafValue = DatumGetKmerP(in->leafDatum);

    if (DatumGetPointer(in->reconstructedValue))
        reconstrValue = (Kmer *) DatumGetPointer(in->reconstructedValue);

    Assert(reconstrValue == NULL ? level == 0 :
           KMER_LEN(reconstrValue) == level);

    /* Reconstruct the full kmer */
    fullLen = level + KMER_LEN(leafValue);
    if (KMER_LEN(leafValue) == 0 && level > 0)
    {
        fullValue = KMER_DATA(reconstrValue);
        out->leafValue = PointerGetDatum(reconstrValue);
    }
    else
    {
        Kmer       *fullKmer = palloc(VARHDRSZ + sizeof(uint8_t) + fullLen);

        SET_VARSIZE(fullKmer, VARHDRSZ + sizeof(uint8_t) + fullLen);
        fullKmer->length = fullLen;
        fullValue = KMER_DATA(fullKmer);
        if (level)
            memcpy(fullValue, KMER_DATA(reconstrValue), level);
        if (KMER_LEN(leafValue) > 0)
            memcpy(fullValue + level, KMER_DATA(leafValue),
                   KMER_LEN(leafValue));
        out->leafValue = PointerGetDatum(fullKmer);
    }

    /* Perform comparisons */
    res = true;
    for (j = 0; j < in->nkeys; j++)
    {
        StrategyNumber strategy = in->scankeys[j].sk_strategy;

        switch (strategy)
        {
            case KMER_EQUAL_STRATEGY:
                {
                    Kmer       *query = DatumGetKmerP(in->scankeys[j].sk_argument);
                    int         queryLen = KMER_LEN(query);

                    if (queryLen != fullLen)
                    {
                        res = false;
                        break;
                    }

                    if (memcmp(fullValue, KMER_DATA(query), fullLen) != 0)
                        res = false;
                }
                break;

            case KMER_PREFIX_STRATEGY:
                {
                    Kmer       *query = DatumGetKmerP(in->scankeys[j].sk_argument);
                    int         queryLen = KMER_LEN(query);

                    /* Check if full kmer starts with query prefix */
                    if (fullLen < queryLen)
                    {
                        res = false;
                        break;
                    }

                    if (memcmp(fullValue, KMER_DATA(query), queryLen) != 0)
                        res = false;
                }
                break;

            case QKMER_CONTAINS_STRATEGY:
                {
                    Qkmer      *query = DatumGetQkmerP(in->scankeys[j].sk_argument);

                    res = qkmer_matches_kmer(QKMER_DATA(query), 
                                             QKMER_LEN(query),
                                             fullValue, 
                                             fullLen);
                }
                break;

            default:
                elog(ERROR, "unrecognized strategy number: %d", strategy);
                res = false;
                break;
        }

        if (!res)
            break;
    }

    PG_RETURN_BOOL(res);
}
