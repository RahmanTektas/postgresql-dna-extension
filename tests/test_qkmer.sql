-- 1. CLEAN STATE (ONLY HERE)
DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;


\echo '========================================='
\echo '===         TEST: QKMER TYPE          ==='
\echo '========================================='


-- Test parsing and output
\echo '\n--- Test 1: Qkmer parsing and output ---'
SELECT 'ACGT'::qkmer AS qk1,
       'ACGT'::qkmer AS qk2,
       'AGCT'::qkmer AS qk3;


-- Test length function
\echo '\n--- Test 2: Qkmer length function ---'
SELECT length('ACGT'::qkmer) AS expected_length_of_qkmer_4;


-- Test IUPAC codes
\echo '\n--- Test 3: IUPAC Codes Parsing (N, R, Y) ---'
SELECT 'ANGT'::qkmer AS wild_N,
       'ARGT'::qkmer AS wild_R,
       'AYGT'::qkmer AS wild_Y;


-- Contains Operator
\echo '\n--- Test 4: Contains Operator (@>) ---'
-- N = Any, R = A or G
SELECT 'ANGT'::qkmer @> 'ACGT'::kmer AS match_N_true,
       'ARGT'::qkmer @> 'AGGT'::kmer AS match_R_true,
       'ARGT'::qkmer @> 'ATGT'::kmer AS match_R_false;


\echo '\n--- Test 5: Multiple ambiguous positions ---'
-- R = A/G, Y = C/T, N = any
SELECT 'RNYT'::qkmer @> 'AGYT'::kmer AS test_same_len_true_1,
       'RNYT'::qkmer @> 'GGCT'::kmer AS test_same_len_true_2,
       'RNYT'::qkmer @> 'CCCT'::kmer AS test_false_wrong_first_base,
       'AANN'::qkmer @> 'AAGT'::kmer AS test_N_any_true,
       'AARY'::qkmer @> 'AAGT'::kmer AS test_RY_true,
       'AARY'::qkmer @> 'AACG'::kmer AS test_RY_false_last_base;


\echo '\n--- Test 6: Full wildcard N-pattern ---'
SELECT 'NNNN'::qkmer @> 'ACGT'::kmer AS full_N_true,
       'NNNN'::qkmer @> 'TTTT'::kmer AS full_N_true2;

-- ======================================================
-- PART 2: DATASET TESTS (FILE BASED)
-- ======================================================

/
\echo '\n--- Test 7: Loading Kmer targets from file ---'
-- We create a table of standard Kmers to test our Qkmer queries against.
DROP TABLE IF EXISTS kmer_targets_dataset;

CREATE TABLE kmer_targets_dataset (
    id bigserial PRIMARY KEY,
    val kmer
);

-- Note: This assumes the file exists. 
-- We load standard Kmers (targets) to filter them later with Qkmers.
COPY kmer_targets_dataset(val)
FROM '/extension/dna_sequence/scriptpy/qkmeroutput.txt'
WITH (FORMAT text);


\echo '\n--- Test 8: Verify Dataset Loading ---'
SELECT count(*) as total_rows_loaded 
FROM kmer_targets_dataset;

SELECT val 
FROM kmer_targets_dataset 
ORDER BY random() 
LIMIT 5;


\echo '\n--- Test 9: Querying Dataset using Wildcard "N" ---'
-- Query: Find all kmers in the dataset that match a pattern with N
-- Example pattern: 'AAAN' (Matches AAAA, AAAC, AAAG, AAAT)
SELECT val 
FROM kmer_targets_dataset 
WHERE 'AAAN'::qkmer @> val
LIMIT 20;


\echo '\n--- Test 10: Querying Dataset using IUPAC Code "R" ---'
-- Query: Find all kmers matching 'AAAR' (Matches AAAA or AAAG only)
SELECT val 
FROM kmer_targets_dataset 
WHERE 'AAAR'::qkmer @> val
LIMIT 20;

\echo '\n--- Test 11: Querying Dataset using IUPAC Code "Y" ---'
SELECT val
FROM kmer_targets_dataset
WHERE 'CYT'::qkmer @> val
LIMIT 20;