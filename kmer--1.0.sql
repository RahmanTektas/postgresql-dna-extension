-- complain if script is sourced in psql, rather than via CREATE EXTENSION
\echo Use "CREATE EXTENSION kmer" to load this file. \quit

/******************************************************************************
 * Input/Output
 ******************************************************************************/

-------------- DNA --------------------------

CREATE OR REPLACE FUNCTION dna_in(cstring)
    RETURNS dna
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_out(dna)
    RETURNS cstring
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/**CREATE OR REPLACE FUNCTION dna_recv(internal)
    RETURNS dna
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;*/

/**CREATE OR REPLACE FUNCTION dna_send(dna)
    RETURNS bytea
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;*/

CREATE TYPE dna (
    internallength = -1, /* -1 représente une taille variable*/
    input = dna_in,
    output = dna_out,
    --receive = dna_recv,
    --send = dna_send,
    alignment = int8,    
    storage = extended,  /*Pour eviter une limite */
);

CREATE CAST (text as dna) WITH FUNCTION dna(text) AS IMPLICIT;
CREATE CAST (dna as text) WITH FUNCTION text(dna);

-------------- Kmer --------------------------

CREATE OR REPLACE FUNCTION kmer_in(cstring)
   RETURNS kmer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_out(kmer)
   RETURNS cstring
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/**CREATE OR REPLACE FUNCTION kmer_recv(internal)
   RETURNS kmer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;*/

/**CREATE OR REPLACE FUNCTION kmer_send(kmer)
   RETURNS bytea
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;*/

CREATE TYPE kmer (
    internallength = 8, 
    input = kmer_in,
    output = kmer_out,
    --receive = kmer_recv,
    --send = kmer_send,
    alignment = int4,    
);

CREATE CAST (text as kmer) WITH FUNCTION kmer(text) AS IMPLICIT;
CREATE CAST (kmer as text) WITH FUNCTION text(kmer);

-------------- Qkmer --------------------------

CREATE OR REPLACE FUNCTION qkmer_in(cstring)
   RETURNS qkmer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_out(qkmer)
   RETURNS cstring
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/**CREATE OR REPLACE FUNCTION qkmer_recv(internal)
   RETURNS qkmer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;*/

/**CREATE OR REPLACE FUNCTION qkmer_send(qkmer)
   RETURNS bytea
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;*/

CREATE TYPE qkmer (
    internallength = 8, 
    input = qkmer_in,
    output = qkmer_out,
   -- receive = qkmer_recv,
   -- send = qkmer_send,
    alignment = int4,    
);

CREATE CAST (text as qkmer) WITH FUNCTION qkmer(text) AS IMPLICIT;
CREATE CAST (qkmer as text) WITH FUNCTION text(qkmer);

/******************************************************************************
 * Constructor
 ******************************************************************************/

-------------- DNA --------------------------
CREATE OR REPLACE FUNCTION dna(text)
  RETURNS dna
  AS 'MODULE_PATHNAME', 'dna_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(dna)
  RETURNS text
  AS 'MODULE_PATHNAME', 'dna_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-------------- Kmer --------------------------
CREATE OR REPLACE FUNCTION kmer(text)
  RETURNS kmer
  AS 'MODULE_PATHNAME', 'kmer_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(kmer)
  RETURNS text
  AS 'MODULE_PATHNAME', 'kmer_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;
  
-------------- Qkmer --------------------------
CREATE OR REPLACE FUNCTION qkmer(text)
  RETURNS qkmer
  AS 'MODULE_PATHNAME', 'qkmer_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(qkmer)
  RETURNS text
  AS 'MODULE_PATHNAME', 'qkmer_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/*****************************************************************************
 * Accessing values
 *****************************************************************************/




/******************************************************************************
 * Operators
 ******************************************************************************/

CREATE FUNCTION  length(dna)
  RETURNS integer
  AS 'MODULE_PATHNAME', 'dna_length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION length(qkmer)
  RETURNS integer
  AS 'MODULE_PATHNAME', 'qkmer_length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION length(kmer)
  RETURNS integer
  AS 'MODULE_PATHNAME', 'kmer_length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION equals(kmer, kmer)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'kmer_eq'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION contains(qkmer, kmer)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'qkmer_contains'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION start_with(kmer, kmer)
  RETURNS boolean  
  AS 'MODULE_PATHNAME', 'kmer_start_with'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION generate_kmers(dna, integer)
  RETURNS SETOF kmer
  AS 'MODULE_PATHNAME', 'generate_kmers'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

  RETURNS SETOF kmer
    AS 'MODULE_PATHNAME', 'generate_kmers'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE

CREATE OPERATOR = (
  LEFTARG = kmer, RIGHTARG = kmer,
  PROCEDURE = equals,
  COMMUTATOR = =, NEGATOR = <>
);

CREATE OPERATOR ^@ (
  LEFTARG = kmer, RIGHTARG = kmer,
  PROCEDURE = start_with
);

CREATE OPERATOR @> (
  LEFTARG = qkmer, RIGHTARG = kmer,
  PROCEDURE = contains
);


