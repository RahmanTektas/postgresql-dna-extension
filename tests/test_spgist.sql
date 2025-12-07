-- 1. CLEAN STATE (ONLY HERE)
DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;


-- 2. SETUP ENVIRONMENT
-- SET enable_seqscan = OFF;


\echo '========================================='
\echo '===        TEST: SP-GIST INDEX        ==='
\echo '========================================='


-- ======================================================
-- PART 1: MANUAL UNIT TESTS (Logic Verification)
-- ======================================================


\echo '\n--- Test 1: Setup and Operator Class Check ---'
-- expecting
--  operator_class_name | access_method | input_type | schema
-- ---------------------+---------------+------------+--------
--  kmer_spgist_ops     | spgist        | kmer       | public

DROP TABLE IF EXISTS kmer_test;
CREATE TABLE kmer_test (val kmer);

SELECT opcname  AS operator_class_name,
       amname   AS access_method,
       opcintype::regtype AS input_type,
       n.nspname AS schema
FROM pg_opclass c
JOIN pg_am a ON c.opcmethod = a.oid
JOIN pg_namespace n ON c.opcnamespace = n.oid
WHERE opcname = 'kmer_spgist_ops';


\echo '\n--- Test 2: Inserting Structured Patterns ---'
INSERT INTO kmer_test VALUES 
    ('AAAAA'), ('AAAAT'), ('AAAAG'), ('AAAAC'),
    ('TGCA'), ('TGCT'), ('TGCG'), ('TGCC'),
    ('CCGG'), ('GGCC'), ('ATAT'), ('TATA');


\echo '\n--- Test 3: Creating SP-GIST Index ---'
CREATE INDEX idx_kmer_spgist
ON kmer_test USING spgist (val public.kmer_spgist_ops);

SET enable_seqscan = OFF;


\echo '\n--- Test 4: Equality Search (=) ---'
EXPLAIN (COSTS OFF)
SELECT val AS matched_kmer
FROM kmer_test
WHERE val = 'TGCA';

SELECT val AS matched_kmer
FROM kmer_test
WHERE val = 'TGCA';


\echo '\n--- Test 5: Prefix Search (^@) ---'
EXPLAIN (COSTS OFF)
SELECT val AS prefix_matched_kmer
FROM kmer_test
WHERE val ^@ 'TGC';

SELECT val AS prefix_matched_kmer
FROM kmer_test
WHERE val ^@ 'TGC';


\echo '--- Test 6: QKmer Pattern Search with "N" (@>) ---'
EXPLAIN (COSTS OFF)
SELECT val AS wildcard_matched_kmer
FROM kmer_test
WHERE 'AAAAN'::qkmer @> val;

SELECT val AS wildcard_matched_kmer
FROM kmer_test
WHERE 'AAAAN'::qkmer @> val;


\echo '--- Test 7: QKmer Pattern Search with "R" (@>) ---'
SELECT val AS wildcard_R_matched_kmer
FROM kmer_test
WHERE 'TGCR'::qkmer @> val;


\echo '\n--- Test 8: Length Mismatch --- : expecting 0 rows'
SELECT val AS length_mismatch_kmer
FROM kmer_test
WHERE 'AAA'::qkmer @> val;



-- ======================================================
-- PART 2: DATASET TESTS 
-- ======================================================

\echo '\n--- Test 12: Create large table for SP-GiST dataset tests ---'
DROP TABLE IF EXISTS spg_kmer_dataset;
CREATE TABLE spg_kmer_dataset(
    id bigserial PRIMARY KEY,
    seq kmer
);

\echo '\n--- Test 13: Loading DNA source file into test_dna ---'
DROP TABLE IF EXISTS test_dna;

CREATE TABLE test_dna(
    id bigserial PRIMARY KEY,
    seq dna
);

COPY test_dna(seq)
FROM '/extension/scriptpy/dna_36.txt'
WITH (FORMAT text);

\echo '\n--- Test 14: Generate kmers (k=5) from DNA file ---'
INSERT INTO spg_kmer_dataset(seq)
SELECT k.kmer AS kmer_value
FROM test_dna d,
     generate_kmers(d.seq, 5) AS k(kmer);


\echo '\n--- Test 15: Create SP-GiST index for dataset ---'
DROP INDEX IF EXISTS spg_kmer_dataset_idx;
CREATE INDEX spg_kmer_dataset_idx
ON spg_kmer_dataset USING spgist (seq);


\echo '\n--- Test 16: Prefix query on dataset ---'
EXPLAIN ANALYZE
SELECT id AS record_id,
       seq AS prefix_matched_kmer
FROM spg_kmer_dataset
WHERE seq ^@ 'ACG'::kmer
LIMIT 10;


\echo '\n--- Test 17: Random equality lookup ---'
EXPLAIN ANALYZE
SELECT id AS record_id,
       seq AS random_matched_kmer
FROM spg_kmer_dataset
WHERE seq = (SELECT seq FROM spg_kmer_dataset ORDER BY random() LIMIT 1);


\echo '\n--- Test 18: qkmer pattern search on dataset ---'
EXPLAIN ANALYZE
SELECT id AS record_id,
       seq AS qkmer_pattern_match
FROM spg_kmer_dataset
WHERE 'RYNNN'::qkmer @> seq
LIMIT 10;


\echo '\n--- Test 19: Count kmers grouped by length ---'
SELECT length(seq) AS kmer_length,
       count(*)    AS count_by_length
FROM spg_kmer_dataset
GROUP BY length(seq)
ORDER BY kmer_length;


\echo '\n--- Test 20: Random sample of dataset ---'
SELECT id  AS record_id,
       seq AS random_sample_kmer
FROM spg_kmer_dataset
ORDER BY random()
LIMIT 5;
