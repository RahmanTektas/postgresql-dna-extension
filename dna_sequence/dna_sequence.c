#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"


#define STB_DS_IMPLEMENTATION
#include "libs/stb_ds.h"

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

/* --- Pointer conversion macros (must come first) --- */
#define DnaPGetDatum(x)   PointerGetDatum(x)
#define DatumGetDnaP(x)   ((Dna *) DatumGetPointer(x))

/* --- Argument and return macros --- */
#define PG_RETURN_DNA_P(x)  return DnaPGetDatum(x)
#define PG_GETARG_DNA_P(n)  DatumGetDnaP(PG_GETARG_DATUM(n))

typedef struct Dna
{
    uint8_t *bases;
    uint8_t length;    /* length of the dna in bases */
} Dna;

/*
 * dna_parse - parse a dna from a string
 *
 * This function takes a pointer to a string and parses it into a Dna structure.
 * If the string is invalid, an error is raised.
 */

static Dna *
dna_parse(const char *str)
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

static char *
dna_to_str(const Dna *dna)
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
Datum
dna_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    Dna *dna = dna_parse(str);   // parse string into dynamic array
    PG_RETURN_DNA_P(dna);
}

PG_FUNCTION_INFO_V1(dna_out);
Datum
dna_out(PG_FUNCTION_ARGS)
{
    Dna *dna = PG_GETARG_DNA_P(0);
    char *str = dna_to_str(dna);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_DNA_P(str);
}

// PG_FUNCTION_INFO_V1(dna_equals);
// Datum
// dna_equals(PG_FUNCTION_ARGS)
// {
//     Dna *a = PG_GETARG_POINTER(0);
//     Dna *b = PG_GETARG_POINTER(1);

//     if (a->length != b->length)
//     {
//         PG_FREE_IF_COPY(a, 0);
//         PG_FREE_IF_COPY(b, 1);
//         PG_RETURN_BOOL(false);
//     }

//     for (uint8_t i = 0; i < a->length; i++)
//     {
//         if (a->bases[i] != b->bases[i])
//         {
//             PG_FREE_IF_COPY(a, 0);
//             PG_FREE_IF_COPY(b, 1);
//             PG_RETURN_BOOL(false);
//         }
//     }

//     PG_FREE_IF_COPY(a, 0);
//     PG_FREE_IF_COPY(b, 1);
//     PG_RETURN_BOOL(true);
// }

PG_FUNCTION_INFO_V1(dna_length);
Datum
dna_length(PG_FUNCTION_ARGS)
{
    Dna *dna = PG_GETARG_DNA_P(0);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_INT32(dna->length);
}

// PG_FUNCTION_INFO_V1(dna_starts_with);
// Datum
// dna_starts_with(PG_FUNCTION_ARGS)
// {
//     Dna *seq = PG_GETARG_POINTER(0);
//     Dna *prefix = PG_GETARG_POINTER(1);

//     if (prefix->length > seq->length)
//         PG_RETURN_BOOL(false);

//     for (uint8_t i = 0; i < prefix->length; i++)
//     {
//         if (seq->bases[i] != prefix->bases[i])
//             PG_RETURN_BOOL(false);
//     }

//     PG_FREE_IF_COPY(seq, 0);
//     PG_FREE_IF_COPY(prefix, 1);
//     PG_RETURN_BOOL(true);
// }



/*  Kmer  */

/* --- Kmer structure --- */
typedef struct Kmer
{
    uint8_t code[32];
    uint8_t  length;    /* length of the kmer in bases */
} Kmer;

/* --- Pointer conversion macros (must come first) --- */
#define KmerPGetDatum(x)   PointerGetDatum(x)
#define DatumGetKmerP(x)   ((Kmer *) DatumGetPointer(x))

/* --- Argument and return macros --- */
#define PG_RETURN_KMER_P(x)  return KmerPGetDatum(x)
#define PG_GETARG_KMER_P(n)  DatumGetKmerP(PG_GETARG_DATUM(n))

// static Kmer *
// kmer_make(uint64_t code, uint8_t length)
// {
//     // allocate zeroed memory for the struct (PostgreSQL’s allocator)
//     Kmer *k = palloc0(sizeof(Kmer));

//     // assign fields
//     k->code = code;
//     k->length = length;

//     // sanity checks
//     if (k->length == 0)
//         ereport(ERROR,
//                 (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
//                  errmsg("k-mer length cannot be zero")));

//     if (k->length > 32)
//         ereport(ERROR,
//                 (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
//                  errmsg("k-mer length exceeds 32 bases")));

//     // all good — return pointer
//     return k;
// }

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



/*  Qkmer  */


/* --- QKmer structure --- */
typedef struct Qkmer
{
    uint8_t code[32];
    uint8_t  length;    /* length of the Qkmer in bases */
} Qkmer;

/* --- Pointer conversion macros (must come first) --- */
#define QkmerPGetDatum(x)   PointerGetDatum(x)
#define DatumGetQkmerP(x)   ((Qkmer *) DatumGetPointer(x))

/* --- Argument and return macros --- */
#define PG_RETURN_QKMER_P(x)  return QkmerPGetDatum(x)
#define PG_GETARG_QKMER_P(n)  DatumGetQkmerP(PG_GETARG_DATUM(n))

static Qkmer *
qkmer_parse(char **str)
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
        switch (s[i])
        {
            case 'A': case 'a': k->code[i] = BASE_A; break;
            case 'C': case 'c': k->code[i] = BASE_C; break;
            case 'G': case 'g': k->code[i] = BASE_G; break;
            case 'T': case 't': k->code[i] = BASE_T; break;
            case 'R': case 'r': k->code[i] = BASE_R; break;
            case 'Y': case 'y': k->code[i] = BASE_Y; break;
            case 'N': case 'n': k->code[i] = BASE_N; break;
            default:
                ereport(ERROR,
                        (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                         errmsg("invalid DNA base in qkmer: \"%c\"", s[i]),
                         errdetail("Only A, C, G, T, R, Y, N are allowed.")));
        }
    }

    return k;
}

static char *
qkmer_to_str(const Qkmer *k)
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
Datum
qkmer_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    PG_RETURN_QKMER_P(qkmer_parse(&str));
}

PG_FUNCTION_INFO_V1(qkmer_out);
Datum
qkmer_out(PG_FUNCTION_ARGS)
{
    Qkmer *k = PG_GETARG_QKMER_P(0);
    char *result = qkmer_to_str(k);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(qkmer_length);
Datum
qkmer_length(PG_FUNCTION_ARGS)
{
    Qkmer *k = PG_GETARG_QKMER_P(0);
    PG_FREE_IF_COPY(k, 0);
    PG_RETURN_INT32(k->length);
}