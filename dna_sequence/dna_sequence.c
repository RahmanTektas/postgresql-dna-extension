#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "dna_sequence.h"
#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"


#define STB_DS_IMPLEMENTATION
#include "libs/stb_ds.h"
PG_MODULE_MAGIC;

//////////////////////////// DNA ////////////////////////////


/*
 * dna_parse - parse a dna from a string
 *
 * This function takes a pointer to a string and parses it into a Dna structure.
 * If the string is invalid, an error is raised.
 */

static Dna * dna_parse(const char *str)
{
    Dna *dna = (Dna *) palloc(sizeof(Dna));
    dna->bases = NULL;
    dna->length = 0;

    size_t len = strlen(str);
    for (size_t i = 0; i < len; i++)
    {
        char c = toupper(str[i]);
        switch (c)
        {
            case 'A': case 'C': case 'G': case 'T':
                arrput(dna->bases, c);  // stb_ds append
                dna->length++;
                break;
            default:
                ereport(ERROR,
                        (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                         errmsg("invalid DNA base: '%c'", str[i]),
                         errdetail("Only A, C, G, T are allowed.")));
        }
    }

    return dna;
}

static char * dna_to_str(const Dna *dna)
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
    char *str = PG_GETARG_CSTRING(0);
    Dna *dna = dna_parse(str);   // parse string into dynamic array
    PG_RETURN_DNA_P(dna);
}

PG_FUNCTION_INFO_V1(dna_out);
Datum dna_out(PG_FUNCTION_ARGS)
    Dna *dna = PG_GETARG_DNA_P(0);
    char *str = dna_to_str(dna);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_CSTRING(str);
}

PG_FUNCTION_INFO_V1(dna_length);
Datum dna_length(PG_FUNCTION_ARGS)
{
    Dna *dna = PG_GETARG_DNA_P(0);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_INT32(dna->length);
}

//////////////////////////// KMER ////////////////////////////


/*
 * kmer_parse - parse a kmer from a string
 *
 * This function takes a pointer to a string and parses it into a Kmer structure.
 * The string is expected to represent a kmer in a specific format (e.g., "ACGT").
 * If the string is invalid, an error is raised.
 */
static Kmer * kmer_parse(char **str)
{
    const char *s = *str;
    int len = strlen(s);

    if (len == 0 || len > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type kmer: \"%s\"", s),
                 errdetail("kmer must be 1–32 bases long")));

    Kmer *k = (Kmer *) palloc0(sizeof(Kmer));
    k->length = len;

    for (int i = 0; i < len; i++)
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
 
static char * kmer_to_str(const Kmer *k)
{
    char *result = palloc(k->length + 1);  /* +1 pour le '\0' */

    for (int i = 0; i < k->length; i++)
    {
        uint8_t mask = k->code[i];

        switch (mask)
        {
            case BASE_A: result[i] = 'A'; break;
            case BASE_C: result[i] = 'C'; break;
            case BASE_G: result[i] = 'G'; break;
            case BASE_T: result[i] = 'T'; break;
            default:     result[i] = '?';  /* au cas où le masque est invalide */
        }
    }

    result[k->length] = '\0';  /* terminaison de chaîne */
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
    Kmer *k = PG_GETARG_KMER_P(0);
    Kmer *j = PG_GETARG_KMER_P(1);

    // If j is longer than k, k cannot start with j
    if (j->length > k->length)
        PG_RETURN_BOOL(false);

    bool result = true;
    for (uint8_t i = 0; i < j->length; i++)
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
        else if (qk->code == BASE_N) continue;                                    

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



}



//////////////////////////// QKMER ////////////////////////////

static Qkmer * qkmer_parse(char **str)
{
    const char *s = *str;
    int len = strlen(s);

    if (len == 0 || len > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type qkmer: \"%s\"", s),
                 errdetail("qkmer must be 1–32 bases long")));

    Qkmer *k = (Qkmer *) palloc0(sizeof(Qkmer));
    k->length = len;

    for (int i = 0; i < len; i++)
    {
        char c = toupper(s[i]);
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

static char * qkmer_to_str(const Qkmer *k)
{
    char *result = palloc(k->length + 1);  /* +1 pour le '\0' */

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
            default:     result[i] = '?';  /* au cas où le masque est invalide */
        }
    }

    result[k->length] = '\0';  /* terminaison de chaîne */
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