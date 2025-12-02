-- 1. CLEAN STATE (ONLY HERE)
DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;


-- 2. SETUP ENVIRONMENT
-- We disable sequential scans to force PostgreSQL to use our SP-GIST index.
-- This is strictly for testing/demonstration purposes on small datasets.
-- SET enable_seqscan = OFF;


\echo '========================================='
\echo '===        TEST: SP-GIST INDEX        ==='
\echo '========================================='


-- ======================================================
-- PART 1: MANUAL UNIT TESTS (Logic Verification)
-- ======================================================


\echo '\n--- Test 1: Setup and Operator Class Check ---'
DROP TABLE IF EXISTS kmer_test;
CREATE TABLE kmer_test (val kmer);

SELECT opcname, amname, opcintype::regtype, n.nspname AS schema
FROM pg_opclass c
JOIN pg_am a ON c.opcmethod = a.oid
JOIN pg_namespace n ON c.opcnamespace = n.oid
WHERE opcname = 'kmer_spgist_ops';


\echo '\n--- Test 2: Inserting Structured Patterns ---'
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


\echo '\n--- Test 3: Creating SP-GIST Index ---'
CREATE INDEX idx_kmer_spgist
ON kmer_test
USING spgist (val public.kmer_spgist_ops);


-- Analyze to update statistics for the query planner
ANALYZE kmer_test;


SET enable_seqscan = OFF;


\echo '\n--- Test 4: Equality Search (=) ---'
-- Should return exactly one row 'TGCA'
SELECT 'Test 4: Equality' as test_name;
EXPLAIN (COSTS OFF) SELECT * FROM kmer_test WHERE val = 'TGCA';
SELECT * FROM kmer_test WHERE val = 'TGCA';


\echo '\n--- Test 5: Prefix Search (^@) ---'
-- Should return 4 rows (TGCA, TGCT, TGCG, TGCC)
-- This tests if the index can match the "inner node" prefix and dump all children.
SELECT 'Test 5: Prefix' as test_name;
EXPLAIN (COSTS OFF) SELECT * FROM kmer_test WHERE val ^@ 'TGC';
SELECT * FROM kmer_test WHERE val ^@ 'TGC';


\echo '--- Test 6: QKmer Pattern Search with "N" (@>) ---'
-- Query: 'AAAAN'
-- Should match: AAAAA, AAAAT, AAAAG, AAAAC
-- Logic: The index must traverse A->A->A->A->[All 4 branches]
SELECT 'Test 6: N Wildcard' as test_name;
EXPLAIN (COSTS OFF) SELECT * FROM kmer_test WHERE 'AAAAN'::qkmer @> val;
SELECT * FROM kmer_test WHERE 'AAAAN'::qkmer @> val;


\echo '--- Test 7: QKmer Pattern Search with "R" (@>) ---'
-- Query: 'TGC R'
-- Should match: TGCA, TGCG
-- Should NOT match: TGCT, TGCC
-- Logic: The index must strictly pick only the A and G branches at the 4th level.
SELECT 'Test 7: R Wildcard' as test_name;
SELECT * FROM kmer_test WHERE 'TGCR'::qkmer @> val ;


\echo '\n--- Test 8: Length Mismatch ---'
-- Query: 'AAA' (Length 3) against data of length 4/5
-- Should return nothing (if strict) or handle gracefully.
SELECT 'Test 8: Length Mismatch' as test_name;
SELECT * FROM kmer_test WHERE 'AAA'::qkmer @> val;


-- ======================================================
-- PART 2: DATASET TESTS 
-- ======================================================

\echo '\n--- Test 12: Create large table for SP-GiST dataset tests ---'
DROP TABLE IF EXISTS spg_kmer_dataset;
CREATE TABLE spg_kmer_dataset(
    id bigserial PRIMARY KEY,
    seq kmer
);


\echo '\n--- Test 13: Generate kmers (k=5) from DNA file ---'
INSERT INTO spg_kmer_dataset(seq)
SELECT k.kmer
FROM test_dna d,
     generate_kmers(d.seq, 5) AS k(kmer);


\echo '\n--- Test 14: Create SP-GiST index for dataset ---'
DROP INDEX IF EXISTS spg_kmer_dataset_idx;
CREATE INDEX spg_kmer_dataset_idx
ON spg_kmer_dataset USING spgist (seq);


\echo '\n--- Test 15: Prefix query on dataset ---'
EXPLAIN ANALYZE
SELECT *
FROM spg_kmer_dataset
WHERE seq ^@ 'ACG'::kmer
LIMIT 10;


\echo '\n--- Test 16: Random equality lookup ---'
EXPLAIN ANALYZE
SELECT *
FROM spg_kmer_dataset
WHERE seq = (SELECT seq FROM spg_kmer_dataset ORDER BY random() LIMIT 1);


\echo '\n--- Test 17: qkmer pattern search on dataset ---'
EXPLAIN ANALYZE
SELECT *
FROM spg_kmer_dataset
WHERE 'RYNNN'::qkmer @> seq
LIMIT 10;


\echo '\n--- Test 18: Count kmers grouped by length ---'
SELECT length(seq) AS len, count(*)
FROM spg_kmer_dataset
GROUP BY length(seq)
ORDER BY len;


\echo '\n--- Test 19: Random sample of dataset ---'
SELECT *
FROM spg_kmer_dataset
ORDER BY random()
LIMIT 5;