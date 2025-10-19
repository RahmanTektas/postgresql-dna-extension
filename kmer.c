#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"

PG_MODULE_MAGIC;

#define EPSILON         1.0E-06

/* --- Kmer structure --- */
typedef struct Kmer
{
    uint64_t code; /* 2 bits per base, up to 32 bases = 64 bits */
    uint8_t  length;    /* length of the kmer in bases */
} Kmer;

/* --- Pointer conversion macros (must come first) --- */
#define KmerPGetDatum(x)   PointerGetDatum(x)
#define DatumGetKmerP(x)   ((Kmer *) DatumGetPointer(x))

/* --- Argument and return macros --- */
#define PG_RETURN_KMER_P(x)  return KmerPGetDatum(x)
#define PG_GETARG_KMER_P(n)  DatumGetKmerP(PG_GETARG_DATUM(n))

/*****************************************************************************/
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

    /* declare variables up-front to avoid mixed-declarations-and-code warning */
    Kmer *k;
    uint64_t bits = 0;
    int i;

    if (len == 0 || len > 32)
        ereport(ERROR,
                (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                 errmsg("invalid input syntax for type kmer: \"%s\"", s),
                 errdetail("kmer must be 1–32 bases long")));

    for (i = 0; i < len; i++)
    {
        bits <<= 2;  /* make room for next base */
        switch (s[i])
        {
            case 'A': case 'a': bits |= 0b00; break;
            case 'C': case 'c': bits |= 0b01; break;
            case 'G': case 'g': bits |= 0b10; break;
            case 'T': case 't': bits |= 0b11; break;
            default:
                ereport(ERROR,
                        (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                         errmsg("invalid DNA base in kmer: \"%c\"", s[i]),
                         errdetail("Only A, C, G, T are allowed.")));
        }
    }

    /* Allocate and fill struct */
    k = (Kmer *) palloc(sizeof(Kmer));
    k->code = bits;
    k->length = len;

    return k;
}

/*
 * kmer_to_str - convert a Kmer structure back to its string representation
 */
static char *
kmer_to_str(const Kmer *k)
{
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    char *result = palloc(k->length + 1);  /* +1 for null terminator */
    uint64_t bits = k->code;
    int i;

    /* Extract bases in reverse (since we encoded left→right) */
    for (i = (int)k->length - 1; i >= 0; i--)
    {
        result[i] = bases[bits & 0b11];  /* decode last 2 bits */
        bits >>= 2;                      /* move to next base */
    }

    result[k->length] = '\0';  /* null-terminate string */
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

PG_FUNCTION_INFO_V1(equals);
Datum
equals(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    Kmer *j = PG_GETARG_KMER_P(1);

    if (k->length != j->length)
            PG_RETURN_BOOL(false);
    PG_FREE_IF_COPY(k, 0);
    PG_FREE_IF_COPY(j, 1);
    PG_RETURN_BOOL(k->code == j->code);
}

PG_FUNCTION_INFO_V1(length);
Datum
length(PG_FUNCTION_ARGS)
{
    Kmer *k = PG_GETARG_KMER_P(0);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_INT32(k->length);
}

PG_FUNCTION_INFO_V1(starts_with);
Datum
starts_with(PG_FUNCTION_ARGS)
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