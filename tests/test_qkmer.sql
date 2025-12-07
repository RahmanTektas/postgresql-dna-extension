-- 1. CLEAN STATE (ONLY HERE)
DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;


\echo '========================================='
\echo '===         TEST: QKMER TYPE          ==='
\echo '========================================='


-- Test parsing and output
\echo '\n--- Test 1: Qkmer parsing and output ---'
SELECT 'ACGT'::qkmer AS qk1,
       'ATGT'::qkmer AS qk2,
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
SELECT 'ANGT'::qkmer @> 'ACGT'::kmer AS match_N_true,   -- true: N covers C
       'ARGT'::qkmer @> 'AGGT'::kmer AS match_R_true,   -- true: R covers G
       'ARGT'::qkmer @> 'ATGT'::kmer AS match_R_false;  -- false: R cannot match T


\echo '\n--- Test 5: Multiple ambiguous positions ---'
SELECT
  'RNYT'::qkmer @> 'AACT'::kmer AS true_test_same_len_1,        -- all positions satisfy constraints
  'RNYT'::qkmer @> 'GGTT'::kmer AS true_test_same_len_2,        -- also valid
  'RNYT'::qkmer @> 'CCCT'::kmer AS false_test_wrong_first_base, -- R does not match C
  'AANN'::qkmer @> 'AAGT'::kmer AS true_test_N_any,             -- N matches G and T
  'AARY'::qkmer @> 'AAGT'::kmer AS true_test_RY,                -- valid match (R=G, Y=T)
  'AARY'::qkmer @> 'AACT'::kmer AS false_test_RY_last_base;     -- Y cannot match C here



\echo '\n--- Test 6: Full wildcard N-pattern ---'
SELECT 'NNNN'::qkmer @> 'ACGT'::kmer AS full_N_true,
       'NNNN'::qkmer @> 'TTTT'::kmer AS full_N_true2;

-- ======================================================
-- PART 2: DATASET TESTS (FILE BASED)
-- ======================================================


\echo '\n--- Test 7: Loading qkmer targets from file ---'
DROP TABLE IF EXISTS qkmer_targets_dataset;

CREATE TABLE qkmer_targets_dataset (
    id  bigserial PRIMARY KEY,
    val qkmer
);

-- Note: This assumes the file exists. 
-- We load standard Kmers (targets) to filter them later with Qkmers.
COPY qkmer_targets_dataset(val)
FROM '/extension/scriptpy/qkmeroutput.txt'
WITH (FORMAT text);


\echo '\n--- Test 8: Verify Dataset Loading ---'
SELECT count(*) as total_rows_loaded 
FROM qkmer_targets_dataset;


\echo '\n--- Test 9: Search for qkmers that include the pattern "AAAN" (Wildcard) : should error ---'
SELECT val AS qkmer_including_AAAN
FROM qkmer_targets_dataset
WHERE val @> 'AAAN'::kmer
LIMIT 5;


\echo '\n--- Test 10: Check whether qkmers cover the kmer "CGTANY" ---'

DROP TABLE IF EXISTS qkmer_with_cgtany;
CREATE TABLE qkmer_with_cgtany (
    id  bigserial PRIMARY KEY,
    val qkmer
);

INSERT INTO qkmer_with_cgtany (val)
VALUES ('CGTANY'::qkmer);

SELECT val AS qkmer_covering_CGTA
FROM qkmer_with_cgtany
WHERE val @> 'CGTA'::kmer;


\echo '\n--- Test 10: Check whether qkmers cover the kmer "CGTA" ---'
SELECT val AS qkmer_covering_CGTA
FROM qkmer_targets_dataset
WHERE val @> 'CGTA'::kmer
--WHERE contains(val,'CGTA'::kmer)
LIMIT 5;


\echo '\n--- Test 11: Check whether qkmers cover the kmer "TCGA" ---'
-- Query: Finds qkmers in the dataset that contain/cover the exact sequence 'TCGA'
SELECT val AS qkmer_covering_TCGA
FROM qkmer_targets_dataset
WHERE val @> 'TCGA'::kmer
LIMIT 5;
