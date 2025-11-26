#include "dna_sequence.h"

//////////////////////////// DNA ////////////////////////////

static Dna *
dna_alloc(uint8 length)
{
    Size size = offsetof(Dna, bases) + length * sizeof(uint8);
    Dna *dna = (Dna *) palloc0(size);

    SET_VARSIZE(dna, size);
    dna->length = length;

    return dna;
}


/*
 * dna_parse - parse a dna from a string
 *
 * This function takes a pointer to a string and parses it into a Dna structure.
 * If the string is invalid, an error is raised.
 */


Dna *dna_parse(const char *str )
{
    size_t len;
    Dna *dna;         
    size_t i;         
    
    
    len = strlen(str);

    if (len == 0 || len > UINT8_MAX)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type dna: \"%s\"", str),
                 errdetail("dna must be 1-255 bases long")));

    for (i = 0; i < len; i++)
    {
        char c = toupper((unsigned char) str[i]);
        switch (c)
        {
            case 'A': case 'C': case 'G': case 'T':
                break;
            default:
                ereport(ERROR,
                        (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                         errmsg("invalid DNA base: '%c'", str[i]),
                         errdetail("Only A, C, G, T are allowed.")));
        }
    }

    dna = dna_alloc((uint8_t) len);

    for (i = 0; i < len; i++)
    {
        dna->bases[i] = (uint8_t) toupper((unsigned char) str[i]);
    }

    return dna;
}


char * dna_to_str(const Dna *dna)
{
    char *result = palloc(dna->length + 1);  // +1 for '\0'

    for (size_t i = 0; i < dna->length; i++)
    {
        // Ensure only canonical bases are returned
        switch (dna->bases[i])
        {
            case 'A': case 'C': case 'G': case 'T':
                result[i] = dna->bases[i];
                break;
            default:
                // fallback: if somehow an invalid base is stored
                result[i] = 'N';
        }
    }

    result[dna->length] = '\0';  // null terminate
    return result;
}

PG_FUNCTION_INFO_V1(dna_in);
Datum dna_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);   /* input function: arg type cstring */
    Dna *dna = dna_parse(str);
    PG_RETURN_DNA_P(dna);
}


PG_FUNCTION_INFO_V1(dna_out);
Datum dna_out(PG_FUNCTION_ARGS)
{
    Dna *dna = PG_GETARG_DNA_P(0);
    char *str = dna_to_str(dna);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_CSTRING(str);   /* output function retourne un cstring */
}

PG_FUNCTION_INFO_V1(dna_length);
Datum dna_length(PG_FUNCTION_ARGS)
{
    Dna *dna = PG_GETARG_DNA_P(0);
    int32 len = dna->length;
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_INT32(len);
}

//////////////////////////// KMER ////////////////////////////

static Kmer *
kmer_alloc(uint8 length)
{
    Size size = offsetof(Kmer, code) + length * sizeof(uint8);
    Kmer *k = (Kmer *) palloc0(size);

    SET_VARSIZE(k, size);
    k->length = length;

    return k;
}


/*
 * kmer_parse - parse a kmer from a string
 *
 * This function takes a pointer to a string and parses it into a Kmer structure.
 * The string is expected to represent a kmer in a specific format (e.g., "ACGT").
 * If the string is invalid, an error is raised.
 */
Kmer * kmer_parse(char **str)
{
    const char *s;
    int len;
    Kmer *k;     
    int i;      

    s = *str;
    len = strlen(s);

    if (len == 0 || len > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type kmer: \"%s\"", s),
                 errdetail("kmer must be 1–32 bases long")));

    k = kmer_alloc((uint8) len);

    for (i = 0; i < len; i++)
    {
        switch (s[i])
        {
            case 'A': case 'a': k->code[i] = BASE_A; break;
            case 'C': case 'c': k->code[i] = BASE_C; break;
            case 'G': case 'g': k->code[i] = BASE_G; break;
            case 'T': case 't': k->code[i] = BASE_T; break;
            default:
                ereport(ERROR,
                        (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                         errmsg("invalid DNA base in kmer: \"%c\"", s[i]),
                         errdetail("Only A, C, G, T are allowed.")));
        }
    }

    return k;
}
 
char * kmer_to_str(const Kmer *k)
{
    char *result = palloc(k->length + 1);

    for (int i = 0; i < k->length; i++)
    {
        uint8_t mask = k->code[i];

        switch (mask)
        {
            case BASE_A: result[i] = 'A'; break;
            case BASE_C: result[i] = 'C'; break;
            case BASE_G: result[i] = 'G'; break;
            case BASE_T: result[i] = 'T'; break;
            default:     result[i] = '?';  
        }
    }

    result[k->length] = '\0';
    return result;
}

PG_FUNCTION_INFO_V1(kmer_in);
Datum kmer_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    PG_RETURN_KMER_P(kmer_parse(&str));
}

PG_FUNCTION_INFO_V1(kmer_out);
Datum kmer_out(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    char *result = kmer_to_str(k);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(kmer_equals);
Datum kmer_equals(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    Kmer *j = PG_GETARG_KMER_P(1);

    // Different lengths → not equal
    if (k->length != j->length)
    {
        PG_FREE_IF_COPY(k, 0);
        PG_FREE_IF_COPY(j, 1);
        PG_RETURN_BOOL(false);
    }

    // Compare each base
    for (uint8_t i = 0; i < k->length; i++)
    {
        if (k->code[i] != j->code[i])
        {
            PG_FREE_IF_COPY(k, 0);
            PG_FREE_IF_COPY(j, 1);
            PG_RETURN_BOOL(false);
        }
    }

    PG_FREE_IF_COPY(k, 0);
    PG_FREE_IF_COPY(j, 1);
    PG_RETURN_BOOL(true);
}

PG_FUNCTION_INFO_V1(kmer_length);
Datum kmer_length(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_INT32(k->length);
}

PG_FUNCTION_INFO_V1(kmer_starts_with);
Datum kmer_starts_with(PG_FUNCTION_ARGS)
{
   Kmer *k;
    Kmer *j;
    bool result;
    uint8_t i;  

    k = PG_GETARG_KMER_P(0);
    j = PG_GETARG_KMER_P(1);

    /* If j is longer than k, k cannot start with j */
    if (j->length > k->length)
    {
        PG_FREE_IF_COPY(k, 0);
        PG_FREE_IF_COPY(j, 1);
        PG_RETURN_BOOL(false);
    }

    result = true;

    for (i = 0; i < j->length; i++)
    {
        // Compare each base using bitwise AND
        // For strict kmer (A/C/G/T only), can just compare equality
        if (k->code[i] != j->code[i])
        {
            result = false;
            break;
        }
    }

    PG_FREE_IF_COPY(k, 0);
    PG_FREE_IF_COPY(j, 1);
    PG_RETURN_BOOL(result);
}

//////////////////////////// QKMER ////////////////////////////

static Qkmer *
qkmer_alloc(uint8 length)
{
    Size size = offsetof(Qkmer, code) + length * sizeof(uint8);
    Qkmer *k = (Qkmer *) palloc0(size);

    SET_VARSIZE(k, size);
    k->length = length;

    return k;
}


PG_FUNCTION_INFO_V1(qkmer_contains);
Datum qkmer_contains(PG_FUNCTION_ARGS)
{
    Qkmer *qk = PG_GETARG_QKMER_P(0);
    Kmer *k = PG_GETARG_KMER_P(1);
    bool result = true;

    if (qk->length != k->length) return false;

    for (size_t i = 0; i < k->length; i++){
        if(k->code[i] == qk->code[i]) continue;

        else if (qk->code[i] == BASE_R && (k->code[i] == BASE_A || k->code[i] == BASE_G)) continue;
        else if (qk->code[i] == BASE_Y && (k->code[i] == BASE_C || k->code[i] == BASE_T)) continue;
        else if (qk->code[i] == BASE_N) continue;                                    

        result = false;
        break;
    }
    PG_FREE_IF_COPY(qk, 0);
    PG_FREE_IF_COPY(k, 1);
    PG_RETURN_BOOL(result);
}

PG_FUNCTION_INFO_V1(generate_kmers);
Datum generate_kmers(PG_FUNCTION_ARGS)
{
    FuncCallContext *funcctx;
    generate_kmers_fctx *fctx;
    Dna *dna;       
    int k;          
    int call_cntr;  
    Kmer *kmer;    
    int i;          

    if (SRF_IS_FIRSTCALL())
    {
        MemoryContext oldcontext;
        funcctx = SRF_FIRSTCALL_INIT();
        oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

        dna = PG_GETARG_DNA_P(0);
        k = PG_GETARG_INT32(1);

        if (k <= 0 || k > dna->length)
            ereport(ERROR, (errmsg("Invalid k")));

        fctx = palloc(sizeof(generate_kmers_fctx));
        fctx->dna_length = dna->length;
        fctx->k = k;
        fctx->num_kmers = dna->length - k + 1;
        fctx->bases = palloc(dna->length);
        memcpy(fctx->bases, dna->bases, dna->length);

        funcctx->user_fctx = fctx;
        funcctx->max_calls = fctx->num_kmers;

        MemoryContextSwitchTo(oldcontext);
    }

    funcctx = SRF_PERCALL_SETUP();
    fctx = funcctx->user_fctx;

    call_cntr = funcctx->call_cntr;

    if (call_cntr < fctx->num_kmers)
    {
        kmer = kmer_alloc((uint8) fctx->k); // +1 ?

        for (i = 0; i < fctx->k; i++)
        {
            char base = toupper(fctx->bases[call_cntr + i]);
            switch (base)
            {
                case 'A': kmer->code[i] = BASE_A; break;
                case 'C': kmer->code[i] = BASE_C; break;
                case 'G': kmer->code[i] = BASE_G; break;
                case 'T': kmer->code[i] = BASE_T; break;
            }
        }

        SRF_RETURN_NEXT(funcctx, PointerGetDatum(kmer));
    }
    else
    {
        SRF_RETURN_DONE(funcctx);
    }
}

//////////////////////////// QKMER ////////////////////////////

Qkmer * qkmer_parse(char **str)
{
    const char *s;
    int len;
    Qkmer *k;  
    int i;     
    char c;    

    s = *str;
    len = strlen(s);

    if (len == 0 || len > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type qkmer: \"%s\"", s),
                 errdetail("qkmer must be 1–32 bases long")));

    k = qkmer_alloc((uint8_t) len);

    for (i = 0; i < len; i++)
    {
        c = toupper(s[i]);
        switch (c)
        {
            case 'A': k->code[i] = BASE_A; break;
            case 'C': k->code[i] = BASE_C; break;
            case 'G': k->code[i] = BASE_G; break;
            case 'T': k->code[i] = BASE_T; break;
            case 'R': k->code[i] = BASE_R; break;
            case 'Y': k->code[i] = BASE_Y; break;
            case 'N': k->code[i] = BASE_N; break;
            default:
                ereport(ERROR,
                        (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                         errmsg("invalid DNA base in qkmer: \"%c\"", s[i]),
                         errdetail("Only A, C, G, T, R, Y, N are allowed.")));
        }
    }

    return k;
}

char * qkmer_to_str(const Qkmer *k)
{
    char *result = palloc(k->length + 1); 

    for (int i = 0; i < k->length; i++)
    {
        uint8_t mask = k->code[i];

        switch (mask)
        {
            case BASE_A: result[i] = 'A'; break;
            case BASE_C: result[i] = 'C'; break;
            case BASE_G: result[i] = 'G'; break;
            case BASE_T: result[i] = 'T'; break;
            case BASE_R: result[i] = 'R'; break;
            case BASE_Y: result[i] = 'Y'; break;
            case BASE_N: result[i] = 'N'; break;
            default:     result[i] = '?'; 
        }
    }

    result[k->length] = '\0';
    return result;
}

PG_FUNCTION_INFO_V1(qkmer_in);
Datum qkmer_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    PG_RETURN_QKMER_P(qkmer_parse(&str));
}

PG_FUNCTION_INFO_V1(qkmer_out);
Datum qkmer_out(PG_FUNCTION_ARGS)
{
    Qkmer *k = PG_GETARG_QKMER_P(0);
    char *result = qkmer_to_str(k);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(qkmer_length);
Datum qkmer_length(PG_FUNCTION_ARGS)
{
    Qkmer *k = PG_GETARG_QKMER_P(0);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_INT32(k->length);
}


/*
 * kmer_hash - hash function for kmer type
 * Required for hash-based operations like GROUP BY
 */

PG_FUNCTION_INFO_V1(kmer_hash);
Datum
kmer_hash(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);

    uint32 hash = hash_any((unsigned char *) k->code, k->length);

    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_UINT32(hash);
}

static int
get_common_prefix_len(const Kmer *a, const Kmer *b)
{
    int i;
    int min_len = (a->length < b->length) ? a->length : b->length;
    
    for (i = 0; i < min_len; i++)
    {
        if (a->code[i] != b->code[i])
            break;
    }
    return i;
}



PG_FUNCTION_INFO_V1(kmer_cmp);
Datum kmer_cmp(PG_FUNCTION_ARGS)
{
    Kmer *a = PG_GETARG_KMER_P(0);
    Kmer *b = PG_GETARG_KMER_P(1);
    int result = 0;
    
    int minlen = (a->length < b->length) ? a->length : b->length;
    result = memcmp(a->code, b->code, minlen);
    
    if (result == 0) {
        if (a->length < b->length) result = -1;
        else if (a->length > b->length) result = 1;
    }
    
    PG_FREE_IF_COPY(a, 0);
    PG_FREE_IF_COPY(b, 1);
    PG_RETURN_INT32(result);
}

// Fix 1: In spg_kmer_config - Use correct type OIDs
PG_FUNCTION_INFO_V1(spg_kmer_config);
Datum spg_kmer_config(PG_FUNCTION_ARGS)
{
    spgConfigIn *cfgin = (spgConfigIn *) PG_GETARG_POINTER(0);
    spgConfigOut *cfg = (spgConfigOut *) PG_GETARG_POINTER(1);

    // Use the input type itself for prefix and leaf
    cfg->prefixType = cfgin->attType;  // Use kmer type for prefix
    cfg->labelType = INT2OID;          // uint8 doesn't have an OID, use INT2 (smallint)
    cfg->leafType = cfgin->attType;    // Use kmer type for leaf
    cfg->canReturnData = true;
    cfg->longValuesOK = false;
    
    PG_RETURN_VOID();
}

// Fix 2: In spg_kmer_choose - Correct struct fields and C90 compliance
PG_FUNCTION_INFO_V1(spg_kmer_choose);
Datum spg_kmer_choose(PG_FUNCTION_ARGS)
{
    spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);
    Kmer *inKmer = DatumGetKmerP(in->datum);
    Kmer *prefixKmer = DatumGetKmerP(in->prefixDatum);
    int commonLen;
    int i;
    Kmer *newPrefix;
    Kmer *postfix;
    Datum *label;

    // Calculate common prefix length
    commonLen = get_common_prefix_len(inKmer, prefixKmer);

    // Case: Need to split the tuple
    if (commonLen < prefixKmer->length)
    {
        out->resultType = spgSplitTuple;
        out->result.splitTuple.prefixHasPrefix = true;
        
        // Create new shorter prefix
        newPrefix = kmer_alloc((uint8_t) commonLen);
        memcpy(newPrefix->code, prefixKmer->code, commonLen);
        out->result.splitTuple.prefixPrefixDatum = PointerGetDatum(newPrefix);
        
        // Create node label - use INT2 (smallint) to store the uint8
        label = (Datum *) palloc(sizeof(Datum));
        *label = Int16GetDatum((int16) prefixKmer->code[commonLen]);
        out->result.splitTuple.childNodeN = 0;  // First child
        
        // Create postfix
        out->result.splitTuple.postfixHasPrefix = true;
        postfix = kmer_alloc((uint8_t) prefixKmer->length - commonLen - 1);
        if (postfix->length > 0)
            memcpy(postfix->code, &prefixKmer->code[commonLen + 1], postfix->length);
        out->result.splitTuple.postfixPrefixDatum = PointerGetDatum(postfix);

        PG_RETURN_VOID();
    }

    // Case: Prefix matches fully, pick a child
    if (inKmer->length > commonLen)
    {
        uint8 nextChar = inKmer->code[commonLen];
        int16 nodeLabel;

        // Look for matching child
        for (i = 0; i < in->nNodes; i++)
        {
            nodeLabel = DatumGetInt16(in->nodeLabels[i]);
            if ((uint8) nodeLabel == nextChar)
            {
                out->resultType = spgMatchNode;
                out->result.matchNode.nodeN = i;
                out->result.matchNode.levelAdd = 1;  // Use levelAdd, not level
                out->result.matchNode.restDatum = in->datum;
                PG_RETURN_VOID();
            }
        }

        // Child not found, add new node
        out->resultType = spgAddNode;
        out->result.addNode.nodeLabel = Int16GetDatum((int16) nextChar);
        out->result.addNode.nodeN = in->nNodes;  // Position for new node
        PG_RETURN_VOID();
    }

    // Exact match - add as new node with terminator
    out->resultType = spgAddNode;
    out->result.addNode.nodeLabel = Int16GetDatum((int16) 0);
    out->result.addNode.nodeN = in->nNodes;
    
    PG_RETURN_VOID();
}

/*
 * SP-GiST 'picksplit' function.
 * Called when a leaf page is full. Creates a new inner tuple.
 */
PG_FUNCTION_INFO_V1(spg_kmer_picksplit);
Datum spg_kmer_picksplit(PG_FUNCTION_ARGS)
{
    spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
    spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);
    Kmer *k0 = DatumGetKmerP(in->datums[0]);
    int commonLen = k0->length;
    int i;
    Kmer *prefixKmer;

    // 1. Find the longest common prefix among ALL datums in this leaf
    for (i = 1; i < in->nTuples; i++)
    {
        Kmer *ki = DatumGetKmerP(in->datums[i]);
        int tmpLen = get_common_prefix_len(k0, ki);
        if (tmpLen < commonLen)
            commonLen = tmpLen;
    }

    // 2. Allocate the new node prefix
    prefixKmer = kmer_alloc((uint8_t) commonLen);

    memcpy(prefixKmer->code, k0->code, commonLen);
    
    out->hasPrefix = true;
    out->prefixDatum = PointerGetDatum(prefixKmer);
    
    // 3. Allocate bucket arrays (we have max 4 bases + maybe terminator, but alloc safe amount)
    // Since we map specific bytes to nodes, we don't know how many distinct next-bytes exist yet.
    out->nNodes = 0;
    out->nodeLabels = (Datum *) palloc(sizeof(Datum) * in->nTuples); 
    out->mapTuplesToNodes = (int *) palloc(sizeof(int) * in->nTuples);
    out->leafTupleDatums = (Datum *) palloc(sizeof(Datum) * in->nTuples);

    // 4. Sort tuples into buckets based on the character *after* the prefix
    for (i = 0; i < in->nTuples; i++)
    {
        Kmer *ki = DatumGetKmerP(in->datums[i]);
        uint8 nextChar = (commonLen < ki->length) ? ki->code[commonLen] : 0;
        int nodeIdx = -1;
        int j;

        // Check if a node for this char already exists
        for (j = 0; j < out->nNodes; j++)
        {
            uint8 *lbl = (uint8 *) DatumGetPointer(out->nodeLabels[j]);
            if (*lbl == nextChar)
            {
                nodeIdx = j;
                break;
            }
        }

        // Create new node label if not found
        if (nodeIdx < 0)
        {
            uint8 *lbl = (uint8 *) palloc(1);
            *lbl = nextChar;
            out->nodeLabels[out->nNodes] = PointerGetDatum(lbl);
            nodeIdx = out->nNodes;
            out->nNodes++;
        }

        out->mapTuplesToNodes[i] = nodeIdx;
        out->leafTupleDatums[i] = in->datums[i];
    }

    PG_RETURN_VOID();
}


// Fix 3: In spg_kmer_inner_consistent - Fix C90 compliance
PG_FUNCTION_INFO_V1(spg_kmer_inner_consistent);
Datum spg_kmer_inner_consistent(PG_FUNCTION_ARGS)
{
    spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
    spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);
    Kmer *prefixKmer = DatumGetKmerP(in->prefixDatum);
    int i, j, p;
    int node_prefix_len = prefixKmer->length;
    int16 nodeLabel;
    uint8 label;
    bool match;
    int label_idx, label_pos, pos, total_len;
    StrategyNumber strategy;
    Kmer *query;
    Qkmer *qk;
    
    out->nNodes = 0;
    out->nodeNumbers = (int *) palloc(sizeof(int) * in->nNodes);
    
    // Check each child branch
    for (i = 0; i < in->nNodes; i++)
    {
        nodeLabel = DatumGetInt16(in->nodeLabels[i]);
        label = (uint8) nodeLabel;
        match = true;

        // Check all query keys
        for (j = 0; j < in->nkeys; j++)
        {
            strategy = in->scankeys[j].sk_strategy;
            
            // Equality (=)
            if (strategy == 1) 
            {
                query = DatumGetKmerP(in->scankeys[j].sk_argument);
                total_len = in->level + node_prefix_len;
                
                // Check prefix match
                for (p = 0; p < node_prefix_len; p++) {
                    if ((in->level + p) >= query->length || 
                        query->code[in->level + p] != prefixKmer->code[p]) {
                        match = false;
                        break;
                    }
                }
                if (!match) break;

                // Check label match
                label_idx = total_len;
                if (label_idx < query->length) {
                    if (query->code[label_idx] != label) match = false;
                } else if (label != 0) {
                    match = false;
                }
            }
            // Prefix (^@)
            else if (strategy == 2) 
            {
                query = DatumGetKmerP(in->scankeys[j].sk_argument);
                
                // Check prefix match
                for (p = 0; p < node_prefix_len; p++) {
                    if ((in->level + p) < query->length && 
                        query->code[in->level + p] != prefixKmer->code[p]) {
                        match = false;
                        break;
                    }
                }
                if (!match) break;

                // Check label
                label_idx = in->level + node_prefix_len;
                if (label_idx < query->length) {
                    if (query->code[label_idx] != label) match = false;
                }
            }
            // Qkmer Containment (@>)
            else if (strategy == 3) 
            {
                qk = DatumGetQkmerP(in->scankeys[j].sk_argument);
                
                // Check node prefix against qkmer
                for (p = 0; p < node_prefix_len; p++) {
                    pos = in->level + p;
                    if (pos >= qk->length) {
                        match = false;
                        break;
                    }
                    
                    // Bitwise check: does qkmer base contain prefix base?
                    if (!(qk->code[pos] & prefixKmer->code[p])) {
                        match = false;
                        break;
                    }
                }
                if (!match) break;

                // Check label
                label_pos = in->level + node_prefix_len;
                if (label_pos >= qk->length) {
                    match = false;
                } else {
                    if (!(qk->code[label_pos] & label)) match = false;
                }
            }
        }

        if (match)
        {
            out->nodeNumbers[out->nNodes] = i;
            out->nNodes++;
        }
    }
    PG_RETURN_VOID();
}




/*
 * SP-GiST 'leaf_consistent' function.
 * Final check on the exact leaf value.
 */
PG_FUNCTION_INFO_V1(spg_kmer_leaf_consistent);
Datum spg_kmer_leaf_consistent(PG_FUNCTION_ARGS)
{
    spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
    spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);
    Kmer *leafKmer = DatumGetKmerP(in->leafDatum);
    bool res = true;
    int j;

    // We must re-check all queries because the tree traversal might be optimistic
    for (j = 0; j < in->nkeys; j++)
    {
        StrategyNumber strategy = in->scankeys[j].sk_strategy;
        
        if (strategy == 1) // =
        {
            Kmer *query = DatumGetKmerP(in->scankeys[j].sk_argument);
            res = DatumGetBool(DirectFunctionCall2(kmer_equals, 
                    PointerGetDatum(leafKmer), PointerGetDatum(query)));
        }
        else if (strategy == 2) // ^@
        {
            Kmer *query = DatumGetKmerP(in->scankeys[j].sk_argument);
            res = DatumGetBool(DirectFunctionCall2(kmer_starts_with, 
                    PointerGetDatum(leafKmer), PointerGetDatum(query)));
        }
        else if (strategy == 3) // @> (qkmer contains kmer)
        {
            Qkmer *qk = DatumGetQkmerP(in->scankeys[j].sk_argument);
            // Note: Argument order for qkmer_contains is (qkmer, kmer)
            res = DatumGetBool(DirectFunctionCall2(qkmer_contains, 
                    PointerGetDatum(qk), PointerGetDatum(leafKmer)));
        }

        if (!res) break;
    }

    out->leafValue = in->leafDatum;
    PG_RETURN_BOOL(res);
}