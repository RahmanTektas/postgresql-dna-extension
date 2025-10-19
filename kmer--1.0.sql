-- complain if script is sourced in psql, rather than via CREATE EXTENSION
\echo Use "CREATE EXTENSION kmer" to load this file. \quit

/******************************************************************************
 * Input/Output
 ******************************************************************************/

-- Type 1: DNA
-- CREATE FUNCTION dna_in(cstring)
--   RETURNS dna
--   AS 'MODULE_PATHNAME', 'dna_in'
--   LANGUAGE C IMMUTABLE STRICT;

-- CREATE FUNCTION dna_out(dna)
--   RETURNS cstring
--   AS 'MODULE_PATHNAME', 'dna_out'
--   LANGUAGE C IMMUTABLE STRICT;

-- CREATE TYPE dna (
--   internallength = VARIABLE,
--   input = dna_in,
--   output = dna_out,
--   alignment = char
-- );

-- Type 2: KMER

CREATE OR REPLACE FUNCTION kmer_in(cstring)
  RETURNS kmer
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_out(kmer)
  RETURNS cstring
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE kmer (
  internallength = 8,    --8 bytes = 64 bits each kmer can store up to 32 nucleotides (2 bits per nucleotide)
  input = kmer_in,
  output = kmer_out,
  alignment = double
);


-- -- Type 3: QKMER
-- CREATE FUNCTION qkmer_in(cstring)
--   RETURNS qkmer
--   AS 'MODULE_PATHNAME', 'qkmer_in'
--   LANGUAGE C IMMUTABLE STRICT;

-- CREATE FUNCTION qkmer_out(qkmer)
--   RETURNS cstring
--   AS 'MODULE_PATHNAME', 'qkmer_out'
--   LANGUAGE C IMMUTABLE STRICT;

-- CREATE TYPE qkmer (
--   internallength = 16,    --16 bytes = 128 bits each qkmer can store up to 32 nucleotides (4 bits per nucleotide)
--   input = qkmer_in,
--   output = qkmer_out,
--   alignment = double
-- );


/*****************************************************************************
 * Accessing values
 *****************************************************************************/
-- CREATE FUNCTION length(dna)
--   RETURNS integer
--   AS 'MODULE_PATHNAME', 'dna_length'
--   LANGUAGE C IMMUTABLE STRICT;

-- CREATE FUNCTION length(kmer)
--   RETURNS integer
--   AS 'MODULE_PATHNAME', 'kmer_length'
--   LANGUAGE C IMMUTABLE STRICT;

-- CREATE FUNCTION length(kmer)
--   RETURNS integer
--   AS 'MODULE_PATHNAME', 'kmer_length'
--   LANGUAGE C IMMUTABLE STRICT;

