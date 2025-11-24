-- complain if script is sourced in psql, rather than via CREATE EXTENSION
\echo Use "CREATE EXTENSION dna_sequence" to load this file. \quit


-------------- DNA --------------------------

CREATE OR REPLACE FUNCTION dna_in(cstring)
    RETURNS dna
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_out(dna)
    RETURNS cstring
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE dna (
    internallength = VARIABLE,
    input = dna_in,
    output = dna_out,
    alignment = int4,       
    storage = extended
);

CREATE OR REPLACE FUNCTION dna(text)
  RETURNS dna
  AS 'MODULE_PATHNAME', 'dna_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(dna)
  RETURNS text
  AS 'MODULE_PATHNAME', 'dna_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

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

CREATE TYPE kmer (
    internallength = 33,          
    input          = kmer_in,
    output         = kmer_out
);

CREATE OR REPLACE FUNCTION kmer(text)
  RETURNS kmer
  AS 'MODULE_PATHNAME', 'kmer_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(kmer)
  RETURNS text
  AS 'MODULE_PATHNAME', 'kmer_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

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

CREATE TYPE qkmer (
    internallength = 8, 
    input = qkmer_in,
    output = qkmer_out,
    alignment = int4  
);

CREATE OR REPLACE FUNCTION qkmer(text)
  RETURNS qkmer
  AS 'MODULE_PATHNAME', 'qkmer_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(qkmer)
  RETURNS text  
  AS 'MODULE_PATHNAME', 'qkmer_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


CREATE CAST (text as qkmer) WITH FUNCTION qkmer(text) AS IMPLICIT;
CREATE CAST (qkmer as text) WITH FUNCTION text(qkmer);

/******************************************************************************
 * Functions
 ******************************************************************************/

CREATE FUNCTION length(dna)
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
  AS 'MODULE_PATHNAME', 'kmer_equals'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION starts_with(kmer, kmer)
  RETURNS boolean  
  AS 'MODULE_PATHNAME', 'kmer_starts_with'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION generate_kmers(dna, integer)
   RETURNS SETOF kmer
   AS 'MODULE_PATHNAME', 'generate_kmers'
   LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION contains(qkmer, kmer)
   RETURNS boolean
   AS 'MODULE_PATHNAME', 'qkmer_contains'
   LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;
   


/******************************************************************************
  * Operators
******************************************************************************/

CREATE OPERATOR = (
  LEFTARG = kmer, RIGHTARG = kmer,
  PROCEDURE = equals,
  COMMUTATOR = =, NEGATOR = <>
);

CREATE OPERATOR ^@ (
  LEFTARG = kmer, RIGHTARG = kmer,
  PROCEDURE = starts_with
);

CREATE OPERATOR @> (
   LEFTARG = qkmer, RIGHTARG = kmer,
   PROCEDURE = contains
 );


/******************************************************************************
 * Hash Support for GROUP BY
******************************************************************************/

CREATE FUNCTION hash(kmer)
  RETURNS integer
  AS 'MODULE_PATHNAME', 'kmer_hash'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- Create hash operator class to enable GROUP BY
CREATE OPERATOR CLASS kmer_hash_ops
    DEFAULT FOR TYPE kmer USING hash AS
    OPERATOR 1 = ,
    FUNCTION 1 hash(kmer);


/******************************************************************************
 * Comparison Support for ORDER BY (BTree)
******************************************************************************/

CREATE FUNCTION kmer_cmp(kmer, kmer)
  RETURNS integer
  AS 'MODULE_PATHNAME', 'kmer_cmp'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- Define comparison operators for btree
CREATE FUNCTION kmer_lt(kmer, kmer)
  RETURNS boolean
  AS 'SELECT kmer_cmp($1, $2) < 0'
  LANGUAGE SQL IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_le(kmer, kmer)
  RETURNS boolean
  AS 'SELECT kmer_cmp($1, $2) <= 0'
  LANGUAGE SQL IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_gt(kmer, kmer)
  RETURNS boolean
  AS 'SELECT kmer_cmp($1, $2) > 0'
  LANGUAGE SQL IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_ge(kmer, kmer)
  RETURNS boolean
  AS 'SELECT kmer_cmp($1, $2) >= 0'
  LANGUAGE SQL IMMUTABLE STRICT PARALLEL SAFE;

-- Create comparison operators
CREATE OPERATOR < (
  LEFTARG = kmer,
  RIGHTARG = kmer,
  PROCEDURE = kmer_lt,
  COMMUTATOR = >,
  NEGATOR = >=,
  RESTRICT = scalarltsel,
  JOIN = scalarltjoinsel
);

CREATE OPERATOR <= (
  LEFTARG = kmer,
  RIGHTARG = kmer,
  PROCEDURE = kmer_le,
  COMMUTATOR = >=,
  NEGATOR = >,
  RESTRICT = scalarlesel,
  JOIN = scalarlejoinsel
);

CREATE OPERATOR > (
  LEFTARG = kmer,
  RIGHTARG = kmer,
  PROCEDURE = kmer_gt,
  COMMUTATOR = <,
  NEGATOR = <=,
  RESTRICT = scalargtsel,
  JOIN = scalargtjoinsel
);

CREATE OPERATOR >= (
  LEFTARG = kmer,
  RIGHTARG = kmer,
  PROCEDURE = kmer_ge,
  COMMUTATOR = <=,
  NEGATOR = <,
  RESTRICT = scalargesel,
  JOIN = scalargejoinsel
);

CREATE FUNCTION kmer_ne(kmer, kmer)
  RETURNS boolean
  AS 'SELECT NOT equals($1, $2)'
  LANGUAGE SQL IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR <> (
  LEFTARG = kmer,
  RIGHTARG = kmer,
  PROCEDURE = kmer_ne,
  COMMUTATOR = <>,
  NEGATOR = =,
  RESTRICT = neqsel,
  JOIN = neqjoinsel
);

/******************************************************************************
 * BTree Operator Class (Required for ORDER BY)
******************************************************************************/

CREATE OPERATOR CLASS kmer_btree_ops
  DEFAULT FOR TYPE kmer USING btree AS
    OPERATOR 1 <,
    OPERATOR 2 <=,
    OPERATOR 3 =,
    OPERATOR 4 >=,
    OPERATOR 5 >,
    FUNCTION 1 kmer_cmp(kmer, kmer);


/******************************************************************************
 * SP-GiST Index Support
******************************************************************************/
-- Register SP-GiST Support Functions
CREATE FUNCTION spg_kmer_config(internal, internal)
  RETURNS void
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION spg_kmer_choose(internal, internal)
  RETURNS void
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION spg_kmer_picksplit(internal, internal)
  RETURNS void
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION spg_kmer_inner_consistent(internal, internal)
  RETURNS void
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION spg_kmer_leaf_consistent(internal, internal)
  RETURNS bool
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT;

-- Define SP-GiST Operator Class
CREATE OPERATOR CLASS kmer_spgist_ops
  FOR TYPE kmer USING spgist AS
    OPERATOR 1  =  (kmer, kmer),
    OPERATOR 2  ^@ (kmer, kmer),
    OPERATOR 3  @> (qkmer, kmer),
    FUNCTION 1  spg_kmer_config(internal, internal),
    FUNCTION 2  spg_kmer_choose(internal, internal),
    FUNCTION 3  spg_kmer_picksplit(internal, internal),
    FUNCTION 4  spg_kmer_inner_consistent(internal, internal),
    FUNCTION 5  spg_kmer_leaf_consistent(internal, internal);