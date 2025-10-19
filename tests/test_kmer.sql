--psql -U postgres -d test_dna -f tests/test_kmer.sql to run the file

-- Drop & recreate the extension to ensure a clean state
DROP EXTENSION IF EXISTS kmer CASCADE;
CREATE EXTENSION kmer;

-- Test parsing and output
SELECT 'ACGT'::kmer AS k1,
       'ACGT'::kmer AS k2,
       'AGCT'::kmer AS k3;

-- Test equality function and operator
SELECT equals('ACGT'::kmer, 'ACGT'::kmer) AS should_be_true;
SELECT equals('ACGT'::kmer, 'AGCT'::kmer) AS should_be_false;
SELECT 'ACGT'::kmer = 'ACGT'::kmer AS operator_true;
SELECT 'ACGT'::kmer = 'AGCT'::kmer AS operator_false;

-- Test length function
SELECT length('ACGT'::kmer) AS length_of_kmer;

-- Test using kmer in a table
DROP TABLE IF EXISTS kmers;
CREATE TABLE kmers (id serial, seq kmer);

INSERT INTO kmers (seq)
VALUES ('ACGT'::kmer), ('AAAA'::kmer), ('AGCT'::kmer);

SELECT * FROM kmers WHERE seq = 'ACGT'::kmer;