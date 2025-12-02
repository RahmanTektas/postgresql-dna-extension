#ifndef DNA_SEQUENCE_H
#define DNA_SEQUENCE_H

#include "postgres.h"
 
#include "access/spgist.h"
#include "access/spgist_private.h"
#include "catalog/pg_type.h"
#include "utils/datum.h"
#include "utils/pg_locale.h"
#include "utils/varlena.h"
#include "common/int.h"
#include "mb/pg_wchar.h"
#include "utils/fmgrprotos.h"
#include "varatt.h"



#include "fmgr.h"
#include <stdint.h> 
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "libpq/pqformat.h"
#include "funcapi.h"  
#include "utils/builtins.h"
#include "access/hash.h"


/* --- Pointer conversion macros --- */
#define DnaPGetDatum(x)   PointerGetDatum(x)
#define DatumGetDnaP(x)   ((Dna *) DatumGetPointer(x))
#define KmerPGetDatum(x)   PointerGetDatum(x)
#define DatumGetKmerP(x)   ((Kmer *) DatumGetPointer(x))
#define QkmerPGetDatum(x)   PointerGetDatum(x)
#define DatumGetQkmerP(x)   ((Qkmer *) DatumGetPointer(x))

/* --- Argument and return macros --- */
#define PG_GETARG_DNA_P(n) \
    ((Dna *) PG_DETOAST_DATUM(PG_GETARG_DATUM(n)))
#define PG_RETURN_DNA_P(x) PG_RETURN_POINTER(x)
#define PG_GETARG_KMER_P(n) \
    ((Kmer *) PG_DETOAST_DATUM(PG_GETARG_DATUM(n)))
#define PG_RETURN_KMER_P(x)  return KmerPGetDatum(x)
#define PG_GETARG_QKMER_P(n) \
    ((Qkmer *) PG_DETOAST_DATUM(PG_GETARG_DATUM(n)))
#define PG_RETURN_QKMER_P(x)  return QkmerPGetDatum(x)

/*  Basic and IUPAC bases  */
#define BASE_A 0x01 // A = 0001
#define BASE_C 0x02 // C = 0010
#define BASE_G 0x04 // G = 0100
#define BASE_T 0x08 // T = 1000

#define BASE_U 0x016
#define BASE_Y (BASE_C | BASE_T)
#define BASE_M (BASE_A | BASE_C)
#define BASE_R (BASE_A | BASE_G)
#define BASE_W (BASE_A | BASE_T)
#define BASE_S (BASE_C | BASE_G)
#define BASE_K (BASE_G | BASE_T)
#define BASE_V (BASE_A | BASE_C | BASE_G)
#define BASE_H (BASE_A | BASE_C | BASE_T)
#define BASE_D (BASE_A | BASE_G | BASE_T)
#define BASE_B (BASE_C | BASE_G | BASE_T)
#define BASE_N (BASE_A | BASE_C | BASE_G | BASE_T)
#define UNKNOWN_SYMBOL '?'

/* --- Dna structure --- */
typedef struct
{
    int32 vl_len_;
    uint8_t length;
    uint8_t bases[FLEXIBLE_ARRAY_MEMBER];
} Dna;

/* --- Kmer structure --- */
typedef struct
{
    int32       vl_len_;        /* varlena header */
    uint8_t       length;         /* 1–255 */
    uint8_t       code[FLEXIBLE_ARRAY_MEMBER];
} Kmer;

/* --- QKmer structure --- */
typedef struct
{
    int32       vl_len_;        /* varlena header */
    uint8_t       length;         /* 1–255 */
    uint8_t       code[FLEXIBLE_ARRAY_MEMBER];
} Qkmer;

typedef struct {
    char *bases;
    int dna_length;
    int k;
    int num_kmers;
} generate_kmers_fctx;

/* Struct for sorting values in picksplit */
typedef struct spgNodePtr
{
    Datum       d;
    int         i;
    int16       c;
} spgNodePtr;

/* Function declarations */
Dna    *dna_parse(const char *str);
char   *dna_to_str(const Dna *dna);
Kmer   *kmer_parse(char **str);
char   *kmer_to_str(const Kmer *k);
Qkmer  *qkmer_parse(char **str);
char   *qkmer_to_str(const Qkmer *k);

#endif
