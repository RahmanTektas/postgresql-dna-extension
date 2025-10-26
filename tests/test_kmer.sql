-- psql -U postgres -d test_dna -f tests/test_kmer.sql to run the file

-- Drop & recreate the extension to ensure a clean state
DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;

\echo '========================================='
\echo '=== KMER TESTS ==========================='
\echo '========================================='

-- Test parsing and output
\echo '\n--- Test 1: Kmer parsing and output ---'
SELECT 'ACGT'::kmer AS k1,
       'ACGT'::kmer AS k2,
       'AGCT'::kmer AS k3;

-- Test equality function and operator
\echo '\n--- Test 2: Kmer equality function ---'
SELECT equals('ACGT'::kmer, 'ACGT'::kmer) AS equality_expected_true;
SELECT equals('ACGT'::kmer, 'AGCT'::kmer) AS equality_expected_false;

\echo '\n--- Test 3: Kmer equality operator ---'
SELECT 'ACGT'::kmer = 'ACGT'::kmer AS equality_operator_true;
SELECT 'ACGT'::kmer = 'AGCT'::kmer AS equality_operator_false;

-- Test length function
\echo '\n--- Test 4: Kmer length ---'
SELECT length('ACGT'::kmer) AS expected_length_of_kmer_4;

-- Tests starts_with(kmer, kmer) function
\echo '\n--- Test 5: Kmer starts_with function ---'
SELECT starts_with('ACGT'::kmer, 'A'::kmer)      AS starts_with_expected_true;
SELECT starts_with('ACGT'::kmer, 'ACGTA'::kmer)  AS starts_with_expected_false;

-- Same tests using operator syntax
\echo '\n--- Test 6: Kmer starts_with operator ---'
SELECT 'ACGT'::kmer ^@ 'A'::kmer      AS starts_with_op_test1_true;
SELECT 'ACGT'::kmer ^@ 'CG'::kmer     AS starts_with_op_test2_false;

-- Test using kmer in a table
\echo '\n--- Test 7: Kmer in table ---'
DROP TABLE IF EXISTS kmers;
CREATE TABLE kmers (id serial, seq kmer);

INSERT INTO kmers (seq)
VALUES ('ACGT'::kmer), ('AAAA'::kmer), ('AGCT'::kmer);

SELECT * FROM kmers WHERE seq = 'ACGT'::kmer;

\echo '\n========================================='
\echo '=== GENERATE_KMERS TESTS ================'
\echo '========================================='

-- Test 1: Basic k=3
\echo '\n--- Test 8: generate_kmers basic (k=3) ---'
SELECT k.kmer
FROM generate_kmers('ACGTACGT', 6) AS k(kmer);

