#include "dna_sequence.h"
#include <stdint.h>


PG_MODULE_MAGIC;

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
    size_t i;

    if (len == 0 || len > UINT8_MAX)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type dna: \"%s\"", str),
                 errdetail("dna must be 1-255 bases long, current length : %li", len)));

    // Validation loop
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
    char *result = palloc(dna->vl_len_ + 1);  // +1 for '\0'

    for (size_t i = 0; i < dna->vl_len_; i++)
    {
        uint8_t mask = char_to_mask((char)dna->bases[i], true);
        if (mask != UNKNOWN_SYMBOL)
        {
             result[i] = (char)dna->bases[i];
        }
        else
        {
             ereport(ERROR,
                    (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                     errmsg("invalid DNA base: '%c'", dna->bases[i]),
                     errdetail("Only A, C, G, T are allowed.")));
        }
    }

    result[dna->vl_len_] = '\0';
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
    size_t i;

    if (len == 0 || len > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type kmer: \"%s\"", s),
                 errdetail("kmer must be 1–32 bases long, current len: %li", len)));

    kmer = kmer_alloc((uint8_t) len);

    for (i = 0; i < len; i++)
    {
        uint8_t mask = char_to_mask(s[i], true);
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
	int len = VARSIZE_ANY_EXHDR(k);
    char *result = palloc(len + 1);

    for (int i = 0; i < len; i++)
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

    if (k_len != j_len)
    {
        PG_FREE_IF_COPY(k, 0);
        PG_FREE_IF_COPY(j, 1);
        PG_RETURN_BOOL(false);
    }

    for (uint8_t i = 0; i < k_len; i++)
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
    PG_RETURN_INT32(VARSIZE_ANY_EXHDR(k));
}

PG_FUNCTION_INFO_V1(kmer_starts_with);
Datum kmer_starts_with(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    Kmer *j = PG_GETARG_KMER_P(1);
	int k_len = VARSIZE_ANY_EXHDR(k);
    int j_len = VARSIZE_ANY_EXHDR(j);
    bool result;
    uint8_t i;
	
    if (j_len > k_len)
    {
        PG_FREE_IF_COPY(k, 0);
        PG_FREE_IF_COPY(j, 1);
        PG_RETURN_BOOL(false);
    }

    result = true;

    for (i = 0; i < j_len; i++)
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

        if (k <= 0 || k > dna->vl_len_)
            ereport(ERROR, (errmsg("Invalid k")));

        fctx = palloc(sizeof(generate_kmers_fctx));
        fctx->dna_length = dna->vl_len_;
        fctx->k = k;
        fctx->num_kmers = dna->vl_len_ - k + 1;
        fctx->bases = palloc(dna->vl_len_);
        memcpy(fctx->bases, dna->bases, dna->vl_len_);

        funcctx->user_fctx = fctx;
        funcctx->max_calls = fctx->num_kmers;

        MemoryContextSwitchTo(oldcontext);
    }

    funcctx = SRF_PERCALL_SETUP();
    fctx = funcctx->user_fctx;
    call_cntr = funcctx->call_cntr;

    if (call_cntr < fctx->num_kmers)
    {
        kmer = kmer_alloc((uint8_t) fctx->k);

        for (i = 0; i < fctx->k; i++)
        {
            // Use helper
            uint8_t mask = char_to_mask((char)fctx->bases[call_cntr + i], true);
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
qkmer_alloc(uint8_t length)
{
    Size size = offsetof(Qkmer, code) + length * sizeof(uint8_t);
    Qkmer *k = (Qkmer *) palloc0(size);
    SET_VARSIZE(k, size);
    return k;
}

Qkmer * qkmer_parse(char **str)
{
    const char *s = *str;
    size_t len = strlen(s);
    Qkmer *k;
    size_t i;

    if (len == 0 || len > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type qkmer: \"%s\"", s),
                 errdetail("qkmer must be 1–32 bases long, current len: %li", len)));

    k = qkmer_alloc((uint8_t) len);

    for (i = 0; i < len; i++)
    {
        // Use helper in non-strict mode (allows all IUPAC)
        uint8_t mask = char_to_mask(s[i], false);
        
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
	int32 len = VARSIZE_ANY_EXHDR(k);
	char *result = palloc(len + 1);

    for (int i = 0; i < len; i++)
    {
        // Use helper
        result[i] = mask_to_char(k->code[i]);
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
    PG_RETURN_INT32(VARSIZE_ANY_EXHDR(k));
}

/*
 * kmer_hash - hash function for kmer type
 */
PG_FUNCTION_INFO_V1(kmer_hash);
Datum
kmer_hash(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    uint32 hash = hash_any((unsigned char *) k->code, VARSIZE_ANY_EXHDR(k));
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
    result = memcmp(a->code, b->code, minlen);
    
    if (result == 0) {
        if (a_len < b_len) result = -1;
        else if (a_len > b_len) result = 1;
    }
    
    PG_FREE_IF_COPY(a, 0);
    PG_FREE_IF_COPY(b, 1);
    PG_RETURN_INT32(result);
}
