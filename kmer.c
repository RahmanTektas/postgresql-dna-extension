#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"

PG_MODULE_MAGIC;

#define EPSILON         1.0E-06

/*  DNA  */

#define BASE_A 0x01 // A = 0001
#define BASE_C 0x02 // C = 0010
#define BASE_G 0x04 // G = 0100
#define BASE_T 0x08 // T = 1000
#define BASE_R (BASE_A | BASE_G) // R = 0001 | 0100
#define BASE_Y (BASE_C | BASE_T) // Y = 0010 | 1000
#define BASE_N (BASE_A | BASE_C | BASE_G | BASE_T) // N = 0001 | 0010 | 0100 | 1000

/* --- Dna structure --- */
// typedef struct Dna
// {
//     uint8_t code;
//     uint8_t  length;    /* length of the dna in bases */
// } Dna;

/*
 * dna_parse - parse a dna from a string
 *
 * This function takes a pointer to a string and parses it into a Dna structure.
 * If the string is invalid, an error is raised.
 */
// static DNA *
// dna_parse(char **str)
// {
//     const char *s = *str;
//     int len = strlen(s);

// }

// PG_FUNCTION_INFO_V1(dna_in);
// Datum
// dna_in(PG_FUNCTION_ARGS)
// {
//     char *str = PG_GETARG_CSTRING(0);
    // PG_RETURN_KMER_P(kmer_parse(&str));
// }


/*  Kmer  */

/* --- Kmer structure --- */
typedef struct Kmer
{
    uint8_t[32] = code;
    // uint64_t code; /* 2 bits per base, up to 32 bases = 64 bits */
    uint8_t  length;    /* length of the kmer in bases */
} Kmer;

/* --- Pointer conversion macros (must come first) --- */
#define KmerPGetDatum(x)   PointerGetDatum(x)
#define DatumGetKmerP(x)   ((Kmer *) DatumGetPointer(x))

/* --- Argument and return macros --- */
#define PG_RETURN_KMER_P(x)  return KmerPGetDatum(x)
#define PG_GETARG_KMER_P(n)  DatumGetKmerP(PG_GETARG_DATUM(n))

static Kmer *
kmer_make(uint64_t code, uint8_t length)
{
    // allocate zeroed memory for the struct (PostgreSQL’s allocator)
    Kmer *k = palloc0(sizeof(Kmer));

    // assign fields
    k->code = code;
    k->length = length;

    // sanity checks
    if (k->length == 0)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("k-mer length cannot be zero")));

    if (k->length > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("k-mer length exceeds 32 bases")));

    // all good — return pointer
    return k;
}

/*
 * kmer_parse - parse a kmer from a string
 *
 * This function takes a pointer to a string and parses it into a Kmer structure.
 * The string is expected to represent a kmer in a specific format (e.g., "ACGT").
 * If the string is invalid, an error is raised.
 */

static Kmer *
kmer_parse(char **str)
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
 
static char *
kmer_to_str(const Kmer *k)
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
Datum
kmer_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    PG_RETURN_KMER_P(kmer_parse(&str));
}

PG_FUNCTION_INFO_V1(kmer_out);
Datum
kmer_out(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    char *result = kmer_to_str(k);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(kmer_equals);
Datum
kmer_equals(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    Kmer *j = PG_GETARG_KMER_P(1);

    if (k->length != j->length)
            PG_RETURN_BOOL(false);
    PG_FREE_IF_COPY(k, 0);
    PG_FREE_IF_COPY(j, 1);
    PG_RETURN_BOOL(k->code == j->code);
}




PG_FUNCTION_INFO_V1(kmer_length);
Datum
kmer_length(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_INT32(k->length);
}

PG_FUNCTION_INFO_V1(kmer_starts_with);
Datum
kmer_starts_with(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    Kmer *j = PG_GETARG_KMER_P(1);
    bool result;
    // If j is longer than k, k cannot start with j
    if (j->length > k->length)
        PG_RETURN_BOOL(false);
    
    // Create a mask for the relevant bits
    uint64_t mask = (1ULL << (2 * j->length)) - 1;
    
    // Shift k's code right to align with j's length
    uint64_t k_prefix = k->code >> (2 * (k->length - j->length));
    
    // Compare the prefixes
    result = ((k_prefix & mask) == j->code);

    PG_FREE_IF_COPY(k, 0);
    PG_FREE_IF_COPY(j, 1);
    PG_RETURN_BOOL(result);
}