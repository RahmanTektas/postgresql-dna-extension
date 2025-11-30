#include "dna_sequence.h"


/*
 * Helper: Convert Char to Bitmask
 * strict: if true, only allows A, C, G, T.
 */
static uint8
char_to_mask(char c, bool strict)
{
    switch (toupper((unsigned char)c))
    {
        // Canonical
        case 'A': return BASE_A;
        case 'C': return BASE_C;
        case 'G': return BASE_G;
        case 'T': return BASE_T;

        // Degenerate (only allowed if not strict)
        case 'M': return strict ? UNKNOWN_SYMBOL : BASE_M;
        case 'R': return strict ? UNKNOWN_SYMBOL : BASE_R;
        case 'W': return strict ? UNKNOWN_SYMBOL : BASE_W;
        case 'S': return strict ? UNKNOWN_SYMBOL : BASE_S;
        case 'Y': return strict ? UNKNOWN_SYMBOL : BASE_Y;
        case 'K': return strict ? UNKNOWN_SYMBOL : BASE_K;
        case 'V': return strict ? UNKNOWN_SYMBOL : BASE_V;
        case 'H': return strict ? UNKNOWN_SYMBOL : BASE_H;
        case 'D': return strict ? UNKNOWN_SYMBOL : BASE_D;
        case 'B': return strict ? UNKNOWN_SYMBOL : BASE_B;
        case 'N': return strict ? UNKNOWN_SYMBOL : BASE_N;

        default: return UNKNOWN_SYMBOL; // Error code
    }
}

/*
 * Helper: Convert Bitmask to Char
 */
static char
mask_to_char(uint8 mask)
{
    switch (mask)
    {
        case BASE_A: return 'A';
        case BASE_C: return 'C';
        case BASE_G: return 'G';
        case BASE_T: return 'T';
        
        // IUPAC
        case BASE_M: return 'M';
        case BASE_R: return 'R';
        case BASE_W: return 'W';
        case BASE_S: return 'S';
        case BASE_Y: return 'Y';
        case BASE_K: return 'K';
        case BASE_V: return 'V';
        case BASE_H: return 'H';
        case BASE_D: return 'D';
        case BASE_B: return 'B';
        case BASE_N: return 'N';
        
        default: return UNKNOWN_SYMBOL;
    }
}

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

Dna *dna_parse(const char *str)
{
    size_t len = strlen(str);
    Dna *dna;
    size_t i;

    if (len == 0 || len > UINT8_MAX)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type dna: \"%s\"", str),
                 errdetail("dna must be 1-255 bases long")));

    // Validation loop
    for (i = 0; i < len; i++)
    {
        uint8 mask = char_to_mask(str[i], true); // true = strict (ACGT only)
        if (mask == UNKNOWN_SYMBOL)
        {
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
        uint8 mask = char_to_mask((char)dna->bases[i], true);
        if (mask != UNKNOWN_SYMBOL)
        {
             result[i] = (char)dna->bases[i]; // It's already stored as char in Dna struct
        }
        else
        {
             ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("invalid DNA base: '%c'", dna->bases[i]),
                     errdetail("Only A, C, G, T are allowed.")));
        }
    }

    result[dna->length] = '\0';
    return result;
}

PG_FUNCTION_INFO_V1(dna_in);
Datum dna_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    Dna *dna = dna_parse(str);
    PG_RETURN_DNA_P(dna);
}


PG_FUNCTION_INFO_V1(dna_out);
Datum dna_out(PG_FUNCTION_ARGS)
{
    Dna *dna = PG_GETARG_DNA_P(0);
    char *str = dna_to_str(dna);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_CSTRING(str);
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

Kmer * kmer_parse(char **str)
{
    const char *s;
    int len;
    Kmer *kmer;
    int i;

    s = *str;
    len = strlen(s);

    if (len == 0 || len > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type kmer: \"%s\"", s),
                 errdetail("kmer must be 1–32 bases long")));

    kmer = kmer_alloc((uint8) len);

    for (i = 0; i < len; i++)
    {
        uint8 mask = char_to_mask(s[i], true); // Strict mode: A, C, G, T only
        if (mask == UNKNOWN_SYMBOL)
        {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("invalid DNA base in kmer: \"%c\"", s[i]),
                     errdetail("Only A, C, G, T are allowed.")));
            kmer->code[i] = UNKNOWN_SYMBOL; 
        }
        else
        {
            kmer->code[i] = mask;
        }
    }

    return kmer;
}

char * kmer_to_str(const Kmer *k)
{
    char *result = palloc(k->length + 1);

    for (int i = 0; i < k->length; i++)
    {
        char c = mask_to_char(k->code[i]);
        if (c == UNKNOWN_SYMBOL)
        {
             ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("invalid internal kmer mask: \"%d\"", k->code[i])));
        }
        result[i] = c;
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

    if (k->length != j->length)
    {
        PG_FREE_IF_COPY(k, 0);
        PG_FREE_IF_COPY(j, 1);
        PG_RETURN_BOOL(false);
    }

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

    if (j->length > k->length)
    {
        PG_FREE_IF_COPY(k, 0);
        PG_FREE_IF_COPY(j, 1);
        PG_RETURN_BOOL(false);
    }

    result = true;

    for (i = 0; i < j->length; i++)
    {
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
        kmer = kmer_alloc((uint8) fctx->k);

        for (i = 0; i < fctx->k; i++)
        {
            // Use helper
            uint8 mask = char_to_mask((char)fctx->bases[call_cntr + i], true);
            if (mask == UNKNOWN_SYMBOL) 
            {
                 // Should not happen for valid Dna type, but safety fallback
                 ereport(ERROR, (errmsg("Invalid base in DNA during kmer generation")));
            }
            kmer->code[i] = mask;
        }

        SRF_RETURN_NEXT(funcctx, PointerGetDatum(kmer));
    }
    else
    {
        SRF_RETURN_DONE(funcctx);
    }
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

    if (qk->length != k->length) {
        PG_FREE_IF_COPY(qk, 0);
        PG_FREE_IF_COPY(k, 1);
        PG_RETURN_BOOL(false);
    } else {
        for (size_t i = 0; i < k->length; i++){
            if(k->code[i] == qk->code[i]) continue;

            else if (qk->code[i] == BASE_R && (k->code[i] == BASE_A || k->code[i] == BASE_G)) continue;
            else if (qk->code[i] == BASE_Y && (k->code[i] == BASE_C || k->code[i] == BASE_T)) continue;
            else if (qk->code[i] == BASE_N) continue;              
            result = false;
            break;
        }
    }
    PG_FREE_IF_COPY(qk, 0);
    PG_FREE_IF_COPY(k, 1);
    PG_RETURN_BOOL(result);
}

Qkmer * qkmer_parse(char **str)
{
    const char *s;
    int len;
    Qkmer *k;
    int i;

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
        // Use helper in non-strict mode (allows all IUPAC)
        uint8 mask = char_to_mask(s[i], false);
        
        if (mask == UNKNOWN_SYMBOL)
        {
            ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("invalid DNA base in qkmer: \"%c\"", s[i]),
                     errdetail("Allowed: A, C, G, T, R, Y, S, W, K, M, B, D, H, V, N.")));
        }
        k->code[i] = mask;
    }

    return k;
}

char * qkmer_to_str(const Qkmer *k)
{
    char *result = palloc(k->length + 1);

    for (int i = 0; i < k->length; i++)
    {
        // Use helper
        result[i] = mask_to_char(k->code[i]);
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

// SP-GiST Choose
PG_FUNCTION_INFO_V1(spg_kmer_choose);
Datum
spg_kmer_choose(PG_FUNCTION_ARGS)
{
    spgChooseIn  *in  = (spgChooseIn *) PG_GETARG_POINTER(0);
    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);

    Kmer *inKmer     = DatumGetKmerP(in->datum);
    Kmer *prefixKmer = DatumGetKmerP(in->prefixDatum);
    int   commonLen;
    int   i;
    Kmer *newPrefix;
    Kmer *postfix;

    commonLen = get_common_prefix_len(inKmer, prefixKmer);

    if (commonLen < prefixKmer->length)
    {
        out->resultType = spgSplitTuple;

        out->result.splitTuple.prefixHasPrefix = true;
        newPrefix = kmer_alloc((uint8) commonLen);
        if (commonLen > 0)
            memcpy(newPrefix->code, prefixKmer->code, commonLen);
        out->result.splitTuple.prefixPrefixDatum = PointerGetDatum(newPrefix);

        out->result.splitTuple.prefixNNodes = 1;
        out->result.splitTuple.prefixNodeLabels = (Datum *) palloc(sizeof(Datum));

        out->result.splitTuple.prefixNodeLabels[0] =
            Int16GetDatum((int16) prefixKmer->code[commonLen]);

        out->result.splitTuple.childNodeN = 0;

        out->result.splitTuple.postfixHasPrefix = true;
        postfix = kmer_alloc((uint8) (prefixKmer->length - commonLen - 1));
        if (postfix->length > 0)
            memcpy(postfix->code,
                   &prefixKmer->code[commonLen + 1],
                   postfix->length);
        out->result.splitTuple.postfixPrefixDatum = PointerGetDatum(postfix);

        PG_RETURN_VOID();
    }

    if (inKmer->length > commonLen)
    {
        uint8 nextChar = inKmer->code[commonLen];
        int16 nodeLabel;

        for (i = 0; i < in->nNodes; i++)
        {
            nodeLabel = DatumGetInt16(in->nodeLabels[i]);
            if ((uint8) nodeLabel == nextChar)
            {
                out->resultType                = spgMatchNode;
                out->result.matchNode.nodeN    = i;
                out->result.matchNode.levelAdd = 1;
                out->result.matchNode.restDatum = in->datum;
                PG_RETURN_VOID();
            }
        }

        out->resultType = spgAddNode;
        out->result.addNode.nodeLabel = Int16GetDatum((int16) nextChar);
        out->result.addNode.nodeN     = in->nNodes;
        PG_RETURN_VOID();
    }

    out->resultType = spgAddNode;
    out->result.addNode.nodeLabel = Int16GetDatum((int16) 0);
    out->result.addNode.nodeN     = in->nNodes;

    PG_RETURN_VOID();
}

// SP-GiST Picksplit
PG_FUNCTION_INFO_V1(spg_kmer_picksplit);
Datum
spg_kmer_picksplit(PG_FUNCTION_ARGS)
{
    spgPickSplitIn  *in  = (spgPickSplitIn *) PG_GETARG_POINTER(0);
    spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);

    Kmer *k0       = DatumGetKmerP(in->datums[0]);
    int   commonLen = k0->length;
    int   i;
    Kmer *prefixKmer;

    for (i = 1; i < in->nTuples; i++)
    {
        Kmer *ki   = DatumGetKmerP(in->datums[i]);
        int   tmpLen = get_common_prefix_len(k0, ki);

        if (tmpLen < commonLen)
            commonLen = tmpLen;
    }

    prefixKmer = kmer_alloc((uint8) commonLen);
    if (commonLen > 0)
        memcpy(prefixKmer->code, k0->code, commonLen);

    out->hasPrefix   = true;
    out->prefixDatum = PointerGetDatum(prefixKmer);

    out->nNodes            = 0;
    out->nodeLabels        = (Datum *) palloc(sizeof(Datum) * in->nTuples);
    out->mapTuplesToNodes = (int *) palloc(sizeof(int) * in->nTuples);
    out->leafTupleDatums  = (Datum *) palloc(sizeof(Datum) * in->nTuples);

    for (i = 0; i < in->nTuples; i++)
    {
        Kmer *ki      = DatumGetKmerP(in->datums[i]);
        uint8 nextChar = (commonLen < ki->length) ? ki->code[commonLen] : 0;
        int   nodeIdx  = -1;
        int   j;

        for (j = 0; j < out->nNodes; j++)
        {
            int16 lbl = DatumGetInt16(out->nodeLabels[j]);
            if ((uint8) lbl == nextChar)
            {
                nodeIdx = j;
                break;
            }
        }

        if (nodeIdx < 0)
        {
            out->nodeLabels[out->nNodes] = Int16GetDatum((int16) nextChar);
            nodeIdx = out->nNodes;
            out->nNodes++;
        }

        out->mapTuplesToNodes[i] = nodeIdx;
        out->leafTupleDatums[i]  = in->datums[i];
    }

    PG_RETURN_VOID();
}

// SP-GiST Inner Consistent
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
    
    for (i = 0; i < in->nNodes; i++)
    {
        nodeLabel = DatumGetInt16(in->nodeLabels[i]);
        label = (uint8) nodeLabel;
        match = true;

        for (j = 0; j < in->nkeys; j++)
        {
            strategy = in->scankeys[j].sk_strategy;
            
            // Equality (=)
            if (strategy == 1) 
            {
                query = DatumGetKmerP(in->scankeys[j].sk_argument);
                total_len = in->level + node_prefix_len;
                
                for (p = 0; p < node_prefix_len; p++) {
                    if ((in->level + p) >= query->length || 
                        query->code[in->level + p] != prefixKmer->code[p]) {
                        match = false;
                        break;
                    }
                }
                if (!match) break;

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
                
                for (p = 0; p < node_prefix_len; p++) {
                    if ((in->level + p) < query->length && 
                        query->code[in->level + p] != prefixKmer->code[p]) {
                        match = false;
                        break;
                    }
                }
                if (!match) break;

                label_idx = in->level + node_prefix_len;
                if (label_idx < query->length) {
                    if (query->code[label_idx] != label) match = false;
                }
            }
            // Qkmer Containment (@>)
            else if (strategy == 3) 
            {
                qk = DatumGetQkmerP(in->scankeys[j].sk_argument);
                
                for (p = 0; p < node_prefix_len; p++) {
                    pos = in->level + p;
                    if (pos >= qk->length) {
                        match = false;
                        break;
                    }
                    
                    // Bitwise check: does qkmer mask contain prefix base mask?
                    if (!(qk->code[pos] & prefixKmer->code[p])) {
                        match = false;
                        break;
                    }
                }
                if (!match) break;

                label_pos = in->level + node_prefix_len;
                if (label_pos >= qk->length) {
                    match = false;
                } else {
                    // Check against the node label (next character in tree)
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

// SP-GiST Leaf Consistent
PG_FUNCTION_INFO_V1(spg_kmer_leaf_consistent);
Datum spg_kmer_leaf_consistent(PG_FUNCTION_ARGS)
{
    spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
    spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);
    Kmer *leafKmer = DatumGetKmerP(in->leafDatum);
    bool res = true;
    int j;

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
        else if (strategy == 3) // @>
        {
            Qkmer *qk = DatumGetQkmerP(in->scankeys[j].sk_argument);
            res = DatumGetBool(DirectFunctionCall2(qkmer_contains, 
                    PointerGetDatum(qk), PointerGetDatum(leafKmer)));
        }

        if (!res) break;
    }

    out->leafValue = in->leafDatum;
    PG_RETURN_BOOL(res);
}
