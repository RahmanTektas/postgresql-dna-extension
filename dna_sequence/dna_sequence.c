#include "dna_sequence.h"
#include <stdint.h>

PG_MODULE_MAGIC;

static inline int
DNA_LEN(Dna *d)
{
    return VARSIZE_ANY_EXHDR(d);
}

static uint8_t
char_to_mask(char c, bool strict)
{
    switch (toupper((unsigned char)c))
    {
        // Default
        case 'A': return BASE_A;
        case 'C': return BASE_C;
        case 'G': return BASE_G;
        case 'T': return BASE_T;
        
        // IUPAC
        case 'U': return strict ? UNKNOWN_SYMBOL : BASE_U;
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

        default: return UNKNOWN_SYMBOL;
    }
}

static char
mask_to_char(uint8_t mask)
{
    switch (mask)
    {
        // Default
        case BASE_A: return 'A';
        case BASE_C: return 'C';
        case BASE_G: return 'G';
        case BASE_T: return 'T';
        
        // IUPAC
        case BASE_U: return 'U';
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
dna_alloc(uint8_t length)
{
    Size size = offsetof(Dna, bases) + length * sizeof(uint8_t);
    Dna *dna = (Dna *) palloc0(size);
    SET_VARSIZE(dna, size);
    return dna;
}

Dna *dna_parse(const char *str)
{
    size_t len = strlen(str);
    Dna *dna;
    // Use VARDATA() for safe write access to the data array
    uint8_t *data; 
    size_t i;

    if (len == 0 || len > UINT8_MAX)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type dna: \"%s\"", str),
                 errdetail("dna must be 1-255 bases long, current length : %li", len)));

    dna = dna_alloc((uint8_t) len);
    data = (uint8_t *) VARDATA(dna); 

    for (i = 0; i < len; i++)
    {
        uint8_t mask = char_to_mask(str[i], true);
        if (mask == UNKNOWN_SYMBOL)
        {
             ereport(ERROR,
                         (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                          errmsg("invalid DNA base: '%c'", str[i]),
                          errdetail("Only A, C, G, T are allowed.")));
        }
        
        data[i] = mask;
    }

    return dna;
}

char * dna_to_str(const Dna *dna)
{
    // Use VARDATA_ANY to safely read the data from a Datum
    uint8_t *data = (uint8_t *) VARDATA_ANY(dna);
	int len = VARSIZE_ANY_EXHDR(dna);
    char *result = palloc(len + 1);  // +1 for '\0'

    for (size_t i = 0; i < len; i++)
    {
        char c = mask_to_char(data[i]);
        result[i] = c;
    }

    result[len] = '\0';
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
    int32 len = VARSIZE_ANY_EXHDR(dna);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_INT32(len);
}

//////////////////////////// KMER ////////////////////////////

static Kmer *
kmer_alloc(uint8_t length)
{
    Size size = offsetof(Kmer, code) + length * sizeof(uint8_t);
    Kmer *k = (Kmer *) palloc0(size);
    SET_VARSIZE(k, size);
    return k;
}

Kmer * kmer_parse(char **str)
{
    const char *s = *str;
    size_t len = strlen(s);
    Kmer *kmer;
    // Use VARDATA() for safe write access
    uint8_t *data; 
    size_t i;

    if (len == 0 || len > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type kmer: \"%s\"", s),
                 errdetail("kmer must be 1–32 bases long, current len: %li", len)));

    kmer = kmer_alloc((uint8_t) len);
    data = (uint8_t *) VARDATA(kmer);

    for (i = 0; i < len; i++)
    {
        uint8_t mask = char_to_mask(s[i], true);
        if (mask == UNKNOWN_SYMBOL)
        {
            ereport(ERROR,
                        (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                         errmsg("invalid DNA base in kmer: \"%c\"", s[i]),
                         errdetail("Only A, C, G, T are allowed.")));
        }
        else
        {
            data[i] = mask;
        }
    }

    return kmer;
}

char * kmer_to_str(const Kmer *k)
{
    // Use VARDATA_ANY to safely read the data
    uint8_t *data = (uint8_t *) VARDATA_ANY(k);
	int len = VARSIZE_ANY_EXHDR(k);
    char *result = palloc(len + 1);

    for (int i = 0; i < len; i++)
    {
        char c = mask_to_char(data[i]);
        if (c == UNKNOWN_SYMBOL)
        {
             ereport(ERROR,
                         (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                          errmsg("invalid internal kmer mask: \"%u\"", data[i])));
        }
        result[i] = c;
    }

    result[len] = '\0';
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
	int k_len = VARSIZE_ANY_EXHDR(k);
	int j_len = VARSIZE_ANY_EXHDR(j);
    bool result = false;

    if (k_len == j_len)
    {
        // Use VARDATA_ANY for safe pointer comparison
        if (memcmp(VARDATA_ANY(k), VARDATA_ANY(j), k_len) == 0)
        {
            result = true;
        }
    }

    PG_FREE_IF_COPY(k, 0);
    PG_FREE_IF_COPY(j, 1);
    PG_RETURN_BOOL(result);
}

PG_FUNCTION_INFO_V1(kmer_length);
Datum kmer_length(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    int32 len = VARSIZE_ANY_EXHDR(k);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_INT32(len);
}

PG_FUNCTION_INFO_V1(kmer_starts_with);
Datum kmer_starts_with(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    Kmer *j = PG_GETARG_KMER_P(1);
	int k_len = VARSIZE_ANY_EXHDR(k);
    int j_len = VARSIZE_ANY_EXHDR(j);
    bool result = false;
	
    if (j_len <= k_len)
    {
        // Use VARDATA_ANY for safe pointer comparison
        if (memcmp(VARDATA_ANY(k), VARDATA_ANY(j), j_len) == 0)
        {
            result = true;
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
	uint8_t *kmer_data;
    int i;

    if (SRF_IS_FIRSTCALL())
    {
        MemoryContext oldcontext;
        funcctx = SRF_FIRSTCALL_INIT();
        oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

        dna = PG_GETARG_DNA_P(0);
        k = PG_GETARG_INT32(1);

        if (k <= 0 || k > 32)
            ereport(ERROR, (errmsg("kmer can not be <= 0 or > 32")));

        if (k > DNA_LEN(dna))
            ereport(ERROR, 
           (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
           errmsg("k cannot be larger than the DNA sequence length")));

        fctx = palloc(sizeof(generate_kmers_fctx));
        fctx->dna_length = DNA_LEN(dna);
        fctx->k = k;
        fctx->num_kmers = fctx->dna_length - k + 1;
        fctx->bases = (uint8_t*) palloc(fctx->dna_length);
        
        // Use VARDATA_ANY to safely copy data from the Dna struct
        memcpy(fctx->bases, VARDATA_ANY(dna), fctx->dna_length);

        funcctx->user_fctx = fctx;
        funcctx->max_calls = fctx->num_kmers;

        MemoryContextSwitchTo(oldcontext);
    }

    funcctx = SRF_PERCALL_SETUP();
    fctx = funcctx->user_fctx;
    call_cntr = funcctx->call_cntr;

    if (call_cntr < fctx->num_kmers)
    {
        kmer = kmer_alloc(fctx->k);
        kmer_data = (uint8_t *) VARDATA(kmer);

        for (i = 0; i < fctx->k; i++)
        {
            uint8_t mask = fctx->bases[call_cntr + i];
            if (mask == UNKNOWN_SYMBOL) 
            {
                ereport(ERROR, (errmsg("Invalid base in DNA during kmer generation")));
            }
            // Write to the safe data pointer
            kmer_data[i] = mask;
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
qkmer_alloc(uint8_t length)
{
    Size size = offsetof(Qkmer, code) + length * sizeof(uint8_t);
    Qkmer *qk = (Qkmer *) palloc0(size);
    SET_VARSIZE(qk, size);
    return qk;
}

Qkmer * qkmer_parse(char **str)
{
    const char *s = *str;
    size_t len = strlen(s);
    Qkmer *qk;
    // Use VARDATA() for safe write access
    uint8_t *data; 
    size_t i;

    if (len == 0 || len > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type qkmer: \"%s\"", s),
                 errdetail("qkmer must be 1–32 bases long, current len: %li", len)));

    qk = qkmer_alloc((uint8_t) len);
    data = (uint8_t *) VARDATA(qk);

    for (i = 0; i < len; i++)
    {
        uint8_t mask = char_to_mask(s[i], false);
        
        if (mask == UNKNOWN_SYMBOL)
        {
            ereport(ERROR,
                        (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                         errmsg("invalid DNA base in qkmer: \"%c\"", s[i]),
                         errdetail("Allowed: A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N.")));
        }
        data[i] = mask;
    }

    return qk;
}

char * qkmer_to_str(const Qkmer *qk)
{
    // Use VARDATA_ANY to safely read the data
    uint8_t *data = (uint8_t *) VARDATA_ANY(qk);
	int32 len = VARSIZE_ANY_EXHDR(qk);
	char *result = palloc(len + 1);

    for (int i = 0; i < len; i++)
    {
        result[i] = mask_to_char(data[i]);
    }

    result[len] = '\0';
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
    Qkmer *qk = PG_GETARG_QKMER_P(0);
    char *result = qkmer_to_str(qk);
    PG_FREE_IF_COPY(qk, 0);
    PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(qkmer_length);
Datum qkmer_length(PG_FUNCTION_ARGS)
{
    Qkmer *qk = PG_GETARG_QKMER_P(0);
    PG_FREE_IF_COPY(qk, 0);
    PG_RETURN_INT32(VARSIZE_ANY_EXHDR(qk));
}

/*
 * kmer_hash - hash function for kmer type
 */
PG_FUNCTION_INFO_V1(kmer_hash);
Datum
kmer_hash(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    // Use VARDATA_ANY to safely pass the pointer to hash_any
    uint32 hash = hash_any((unsigned char *) VARDATA_ANY(k), VARSIZE_ANY_EXHDR(k));
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_UINT32(hash);
}

PG_FUNCTION_INFO_V1(kmer_cmp);
Datum kmer_cmp(PG_FUNCTION_ARGS)
{
    Kmer *a = PG_GETARG_KMER_P(0);
    Kmer *b = PG_GETARG_KMER_P(1);
	int a_len = VARSIZE_ANY_EXHDR(a);
    int b_len = VARSIZE_ANY_EXHDR(b);
	
    int result = 0;
    
    int minlen = (a_len < b_len) ? a_len : b_len;
    // Use VARDATA_ANY for safe data comparison
    result = memcmp(VARDATA_ANY(a), VARDATA_ANY(b), minlen);
    
    if (result == 0) {
        if (a_len < b_len) result = -1;
        else if (a_len > b_len) result = 1;
    }
    
    PG_FREE_IF_COPY(a, 0);
    PG_FREE_IF_COPY(b, 1);
    PG_RETURN_INT32(result);
}
