-- Drop & recreate the extension to ensure a clean state
DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;

-- Test parsing and output
SELECT 'ACGT'::qkmer AS qk1,
       'ACGT'::qkmer AS qk2,
       'AGCT'::qkmer AS qk3;

-- Test length function
SELECT length('ACGT'::qkmer) AS expected_length_of_qkmer_4;

-- Test using qkmer in a table
DROP TABLE IF EXISTS qkmers;
CREATE TABLE qkmers (id serial, seq qkmer);