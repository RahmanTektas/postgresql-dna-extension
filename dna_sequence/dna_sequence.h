#ifndef DNA_SEQUENCE_H
#define DNA_SEQUENCE_H

#include <stdint.h> 
#include <math.h>
#include <float.h>
#include <stdlib.h>

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
#define PG_RETURN_DNA_P(x)  return DnaPGetDatum(x)
#define PG_GETARG_DNA_P(n)  DatumGetDnaP(PG_GETARG_DATUM(n))


/* --- Pointer conversion macros (must come first) --- */
#define QkmerPGetDatum(x)   PointerGetDatum(x)
#define DatumGetQkmerP(x)   ((Qkmer *) DatumGetPointer(x))

/* --- Argument and return macros --- */
#define PG_RETURN_QKMER_P(x)  return QkmerPGetDatum(x)
#define PG_GETARG_QKMER_P(n)  DatumGetQkmerP(PG_GETARG_DATUM(n))

/* --- Dna structure --- */
typedef struct Dna
{
    uint8_t *bases;
    uint8_t length;    /* length of the dna in bases */
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

// Prototypes des fonctions internes
Dna    *dna_parse(const char *str);
char   *dna_to_str(const Dna *dna);
Kmer   *kmer_parse(char **str);
char   *kmer_to_str(const Kmer *k);
Qkmer  *qkmer_parse(char **str);
char   *qkmer_to_str(const Qkmer *k);

#endif
