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

CREATE OR REPLACE FUNCTION kmer_out(kmer)
  RETURNS cstring
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE kmer (
  internallength = 8,    --8 bytes = 64 bits each kmer can store up to 32 nucleotides (2 bits per nucleotide)
  input = kmer_in,
  output = kmer_out,
  alignment = double
);

CREATE FUNCTION equals(kmer, kmer)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'equals'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR = (
  LEFTARG = kmer,
  RIGHTARG = kmer,
  PROCEDURE = equals,
  COMMUTATOR = =,     --says that a = b is the same as b = a
  NEGATOR = <>        --says that the negation of a = b is a <> b
);

CREATE OR REPLACE FUNCTION length(kmer)
  RETURNS integer
  AS 'MODULE_PATHNAME', 'length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION starts_with(kmer, kmer)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'starts_with'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

 
CREATE OPERATOR ^@ (
  LEFTARG = kmer,
  RIGHTARG = kmer,
  PROCEDURE = starts_with
);

-- CREATE FUNCTION length(kmer)
--   RETURNS integer
--   AS 'MODULE_PATHNAME', 'kmer_length'
--   LANGUAGE C IMMUTABLE STRICT;



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

