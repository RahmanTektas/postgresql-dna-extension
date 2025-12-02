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
SELECT 'ANGT'::qkmer @> 'ACGT'::kmer AS match_N_true,   -- true: N covers C
       'ARGT'::qkmer @> 'AGGT'::kmer AS match_R_true,   -- true: R covers G
       'ARGT'::qkmer @> 'ATGT'::kmer AS match_R_false;  -- false: R cannot match T


\echo '\n--- Test 5: Multiple ambiguous positions ---'
SELECT
  'RNYT'::qkmer @> 'AACT'::kmer AS test_same_len_true_1,        -- all positions satisfy constraints
  'RNYT'::qkmer @> 'GGTT'::kmer AS test_same_len_true_2,        -- also valid
  'RNYT'::qkmer @> 'CCCT'::kmer AS test_false_wrong_first_base, -- R does not match C
  'AANN'::qkmer @> 'AAGT'::kmer AS test_N_any_true,             -- N matches G and T
  'AARY'::qkmer @> 'AAGT'::kmer AS test_RY_true,                -- valid match (R=G, Y=T)
  'AARY'::qkmer @> 'AACT'::kmer AS test_RY_false_last_base;     -- Y cannot match C here



\echo '\n--- Test 6: Full wildcard N-pattern ---'
SELECT 'NNNN'::qkmer @> 'ACGT'::kmer AS full_N_true,
       'NNNN'::qkmer @> 'TTTT'::kmer AS full_N_true2;

-- ======================================================
-- PART 2: DATASET TESTS (FILE BASED)
-- ======================================================


\echo '\n--- Test 7: Loading qkmer targets from file ---'
DROP TABLE IF EXISTS qkmer_targets_dataset;

CREATE TABLE kmer_targets_dataset (
    id  bigserial PRIMARY KEY,
    val kmer
);

-- Note: This assumes the file exists. 
-- We load standard Kmers (targets) to filter them later with Qkmers.
COPY qkmer_targets_dataset(val)
FROM '/extension/scriptpy/qkmeroutput.txt'
WITH (FORMAT text);


\echo '\n--- Test 8: Verify Dataset Loading ---'
SELECT count(*) as total_rows_loaded 
FROM kmer_targets_dataset;


\echo '\n--- Test 9: Recherche de qkmers incluant le motif "AAAN" (Wildcard) ---'
SELECT val
FROM qkmer_targets_dataset
WHERE val @> 'AAAN'::kmer
LIMIT 20;

\echo '\n--- Test 10: Vérification de couverture du kmer "CGTA" ---'
SELECT val 
FROM qkmer_targets_dataset 
WHERE val @> 'CGTA'::kmer
LIMIT 20;


\echo '\n--- Test 11: Vérification de couverture du kmer "TCGA" ---'
-- Query: Trouve les qkmers du dataset qui contiennent/couvrent la séquence exacte 'TCGA'
SELECT val
FROM qkmer_targets_dataset
WHERE val @> 'TCGA'::kmer
LIMIT 20;