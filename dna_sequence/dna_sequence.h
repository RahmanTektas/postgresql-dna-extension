#ifndef DNA_SEQUENCE_H
#define DNA_SEQUENCE_H

/* PostgreSQL headers MUST come first */
#include "postgres.h"
#include "fmgr.h"
#include "utils/varlena.h"
#include <stdint.h> 
#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include "funcapi.h"  
#include "utils/builtins.h"

PG_MODULE_MAGIC;

/* --- Pointer conversion macros (must come first) --- */
#define KmerPGetDatum(x)   PointerGetDatum(x)
#define DatumGetKmerP(x)   ((Kmer *) DatumGetPointer(x))

/* --- Argument and return macros --- */
#define PG_RETURN_KMER_P(x)  return KmerPGetDatum(x)
#define PG_GETARG_KMER_P(n)  DatumGetKmerP(PG_GETARG_DATUM(n))


/*  DNA  */
#define BASE_A 0x01 // A = 0001
#define BASE_C 0x02 // C = 0010
#define BASE_G 0x04 // G = 0100
#define BASE_T 0x08 // T = 1000
#define BASE_R (BASE_A | BASE_G) // R = 0001 | 0100
#define BASE_Y (BASE_C | BASE_T) // Y = 0010 | 1000
#define BASE_N (BASE_A | BASE_C | BASE_G | BASE_T) // N = 0001 | 0010 | 0100 | 1000


/* --- Pointer conversion macros (must come first) --- */
#define DnaPGetDatum(x)   PointerGetDatum(x)
#define DatumGetDnaP(x)   ((Dna *) DatumGetPointer(x))

/* --- Argument and return macros --- */
#define PG_GETARG_DNA_P(n) ((Dna *) PG_DETOAST_DATUM(PG_GETARG_DATUM(n)))
#define PG_RETURN_DNA_P(x) PG_RETURN_POINTER(x)


/* --- Pointer conversion macros (must come first) --- */
#define QkmerPGetDatum(x)   PointerGetDatum(x)
#define DatumGetQkmerP(x)   ((Qkmer *) DatumGetPointer(x))

/* --- Argument and return macros --- */
#define PG_RETURN_QKMER_P(x)  return QkmerPGetDatum(x)
#define PG_GETARG_QKMER_P(n)  DatumGetQkmerP(PG_GETARG_DATUM(n))

/* --- Dna structure --- */
typedef struct
{
    int32 vl_len_;
    uint8_t length;
    //uint8_t pad[3];  /* padding to keep bases[] 4-byte aligned */
    uint8_t bases[FLEXIBLE_ARRAY_MEMBER];
} Dna;

/* --- Kmer structure --- */
typedef struct Kmer
{
    uint8_t code[32];
    uint8_t  length;    /* length of the kmer in bases */
} Kmer;

/* --- QKmer structure --- */
typedef struct Qkmer
{
    uint8_t code[32];
    uint8_t  length;    /* length of the Qkmer in bases */
} Qkmer;

typedef struct {
    char *bases;
    int dna_length;
    int k;
    int num_kmers;
} generate_kmers_fctx;

/* Function declarations */
Dna    *dna_parse(const char *str);
char   *dna_to_str(const Dna *dna);
Kmer   *kmer_parse(char **str);
char   *kmer_to_str(const Kmer *k);
Qkmer  *qkmer_parse(char **str);
char   *qkmer_to_str(const Qkmer *k);

#endif
