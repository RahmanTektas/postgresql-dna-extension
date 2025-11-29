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


-- ======================================================
-- PART 2: DATASET TESTS (FILE BASED)
-- ======================================================

/**
\echo '\n--- Test 5: Loading Kmer targets from file ---'
-- We create a table of standard Kmers to test our Qkmer queries against.
DROP TABLE IF EXISTS kmer_targets_dataset;

CREATE TABLE kmer_targets_dataset (
    id bigserial PRIMARY KEY,
    val kmer
);

-- Note: This assumes the file exists. 
-- We load standard Kmers (targets) to filter them later with Qkmers.
COPY kmer_targets_dataset(val)
FROM '/extension/dna_sequence/scriptpy/qkmers_targets.txt'
WITH (FORMAT text);


\echo '\n--- Test 6: Verify Dataset Loading ---'
SELECT count(*) as total_rows_loaded 
FROM kmer_targets_dataset;

SELECT val 
FROM kmer_targets_dataset 
ORDER BY random() 
LIMIT 5;


\echo '\n--- Test 7: Querying Dataset using Wildcard "N" ---'
-- Query: Find all kmers in the dataset that match a pattern with N
-- Example pattern: 'AAAN' (Matches AAAA, AAAC, AAAG, AAAT)
SELECT val 
FROM kmer_targets_dataset 
WHERE 'AAAN'::qkmer @> val
LIMIT 20;


\echo '\n--- Test 8: Querying Dataset using IUPAC Code "R" ---'
-- Query: Find all kmers matching 'AAAR' (Matches AAAA or AAAG only)
SELECT val 
FROM kmer_targets_dataset 
WHERE 'AAAR'::qkmer @> val
LIMIT 20;

**/