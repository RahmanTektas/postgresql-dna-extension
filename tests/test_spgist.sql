DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;

-- 1. Setup Environment
DROP TABLE IF EXISTS kmer_test;
CREATE TABLE kmer_test (val kmer);

SELECT opcname, amname, opcintype::regtype, n.nspname AS schema
FROM pg_opclass c
JOIN pg_am a ON c.opcmethod = a.oid
JOIN pg_namespace n ON c.opcnamespace = n.oid
WHERE opcname = 'kmer_spgist_ops';

-- 2. Insert Data
-- We insert specific patterns to test the Trie structure
-- Group 1: Poly-A variations (tests depth)
INSERT INTO kmer_test VALUES 
    ('AAAAA'), ('AAAAT'), ('AAAAG'), ('AAAAC');

-- Group 2: Prefix splitting
-- These share 'TGC' but diverge at the 4th base
INSERT INTO kmer_test VALUES 
    ('TGCA'), ('TGCT'), ('TGCG'), ('TGCC');

-- Group 3: Random noise to fill pages
INSERT INTO kmer_test VALUES 
    ('CCGG'), ('GGCC'), ('ATAT'), ('TATA');

-- 3. Create the SP-GiST Index
CREATE INDEX idx_kmer_spgist
ON kmer_test
USING spgist (val public.kmer_spgist_ops);

-- 4. ANALYZE to update statistics
ANALYZE kmer_test;

-- 5. FORCE INDEX USAGE
-- This is crucial. Without this, PG might just scan the table because it's small.
SET enable_seqscan = OFF;

---------------------------------------------------------------------
-- TEST 1: Equality (=)
-- Should return exactly one row 'TGCA'
---------------------------------------------------------------------
SELECT 'Test 1: Equality' as test_name;
EXPLAIN (COSTS OFF) SELECT * FROM kmer_test WHERE val = 'TGCA';
SELECT * FROM kmer_test WHERE val = 'TGCA';

---------------------------------------------------------------------
-- TEST 2: Prefix Search (^@)
-- Should return 4 rows (TGCA, TGCT, TGCG, TGCC)
-- This tests if the index can match the "inner node" prefix and dump all children.
---------------------------------------------------------------------
SELECT 'Test 2: Prefix' as test_name;
EXPLAIN (COSTS OFF) SELECT * FROM kmer_test WHERE val ^@ 'TGC';
SELECT * FROM kmer_test WHERE val ^@ 'TGC';

---------------------------------------------------------------------
-- TEST 3: Wildcard Containment (@>) - The "N" Wildcard
-- Query: 'AAAAN'
-- Should match: AAAAA, AAAAT, AAAAG, AAAAC
-- Logic: The index must traverse A->A->A->A->[All 4 branches]
---------------------------------------------------------------------
SELECT 'Test 3: N Wildcard' as test_name;
EXPLAIN (COSTS OFF) SELECT * FROM kmer_test WHERE 'AAAAN'::qkmer @> val;
SELECT * FROM kmer_test WHERE 'AAAAN'::qkmer @> val;

---------------------------------------------------------------------
-- TEST 4: Wildcard Containment (@>) - The "R" Wildcard (Purine: A or G)
-- Query: 'TGC R'
-- Should match: TGCA, TGCG
-- Should NOT match: TGCT, TGCC
-- Logic: The index must strictly pick only the A and G branches at the 4th level.
---------------------------------------------------------------------
SELECT 'Test 4: R Wildcard' as test_name;
SELECT * FROM kmer_test WHERE 'TGCR'::qkmer @> val ;

---------------------------------------------------------------------
-- TEST 5: Mismatch length
-- Query: 'AAA' (Length 3) against data of length 4/5
-- Should return nothing (if strict) or handle gracefully.
---------------------------------------------------------------------
SELECT 'Test 5: Length Mismatch' as test_name;
EXPLAIN ANALYZE SELECT * FROM kmer_test WHERE 'AAA'::qkmer @> val;



SELECT * 
FROM pg_class
WHERE relname = 'idx_kmer_spgist';

SELECT * 
FROM spgist_page_items(get_raw_page('idx_kmer_spgist', 0));


WITH pages AS (
  SELECT
    blkno,
    (spgist_page_items(get_raw_page('idx_kmer_spgist', blkno))).*
  FROM generate_series(0,
        (SELECT relpages - 1
         FROM pg_class
         WHERE relname = 'idx_kmer_spgist')) AS blkno
)
SELECT
  level,
  blkno,
  tupleoffset,
  -- prefix is stored as the datum; cast it to kmer for readability
  prefix::kmer AS node_prefix,
  nodeLabel::smallint AS child_label
FROM pages
WHERE type = 'i'      -- inner tuples only
ORDER BY level, blkno, tupleoffset;