#include "dna_sequence.h"

//////////////////////////// DNA ////////////////////////////


/*
 * dna_parse - parse a dna from a string
 *
 * This function takes a pointer to a string and parses it into a Dna structure.
 * If the string is invalid, an error is raised.
 */


Dna *dna_parse(const char *str)
{
    size_t len;
    Size size;        
    Dna *dna;         
    size_t i;         
    
    
    len = strlen(str);

    if (len > UINT8_MAX)
        ereport(ERROR,
                (errcode(ERRCODE_STRING_DATA_RIGHT_TRUNCATION),
                 errmsg("DNA sequence too long (max %d bases)", UINT8_MAX)));

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

    size = offsetof(Dna, bases) + len * sizeof(uint8_t);

    dna = (Dna *) palloc(size);
    SET_VARSIZE(dna, size);

    dna->length = (uint8_t) len;

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

PG_FUNCTION_INFO_V1(dna_cast_from_text);
Datum dna_cast_from_text(PG_FUNCTION_ARGS)
{
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = text_to_cstring(txt);

    Dna *dna = dna_parse(str);

    pfree(str);
    PG_FREE_IF_COPY(txt, 0);
    PG_RETURN_DNA_P(dna);
}

PG_FUNCTION_INFO_V1(dna_cast_to_text);
Datum
dna_cast_to_text(PG_FUNCTION_ARGS)
{
    Dna *dna;
    text *result;
    char *str;
    int i; 

    dna = PG_GETARG_DNA_P(0);

    if (dna == NULL || dna->length == 0)
        PG_RETURN_TEXT_P(cstring_to_text(""));

    /* Allocate buffer for string (length + 1 for null terminator) */
    str = (char *) palloc(dna->length + 1);
    
    /* Copy bases to string */
    for (i = 0; i < dna->length; i++)
    {
        str[i] = dna->bases[i];
    }
    str[dna->length] = '\0';
    
    result = cstring_to_text(str);
    pfree(str);
    
    PG_RETURN_TEXT_P(result);
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

    k = (Kmer *) palloc0(sizeof(Kmer));
    k->length = len;

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

PG_FUNCTION_INFO_V1(kmer_cast_from_text);
Datum kmer_cast_from_text(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    Kmer *k = kmer_parse(&str);
    PG_RETURN_KMER_P(k);
}

PG_FUNCTION_INFO_V1(kmer_cast_to_text);
Datum kmer_cast_to_text(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    char *str = kmer_to_str(k);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_CSTRING(str);
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
        // Allouer et remplir un Kmer
        kmer = palloc0(sizeof(Kmer)+1);
        kmer->length = fctx->k;

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


PG_FUNCTION_INFO_V1(qkmer_cast_from_text);
Datum qkmer_cast_from_text(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    Qkmer *qk = qkmer_parse(&str);
    PG_RETURN_QKMER_P(qk);
}

PG_FUNCTION_INFO_V1(qkmer_cast_to_text);
Datum qkmer_cast_to_text(PG_FUNCTION_ARGS)
{
    Qkmer *qk = PG_GETARG_QKMER_P(0);
    char *str = qkmer_to_str(qk);
    PG_FREE_IF_COPY(qk, 0);
    PG_RETURN_CSTRING(str);
}

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

    k = (Qkmer *) palloc0(sizeof(Qkmer));
    k->length = len;

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
Datum kmer_hash(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    uint32 hash = 0;
    
    // Hash combines length and all bases
    hash = (uint32) k->length;
    
    for (int i = 0; i < k->length; i++)
    {
        // hash = hash * 33 + code[i] (djb2 algorithm variant)
        hash = ((hash << 5) + hash) + (uint32) k->code[i];
    }
    
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_INT32(hash);
}
