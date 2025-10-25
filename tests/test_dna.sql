-- Run this file using:
--   psql -U postgres -d test_dna -f tests/test_dna.sql

-- Clean state
DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;


-- Basic parsing and output tests
SELECT 'ACGT'::dna AS dna1,
       'ACGN'::dna AS dna2,
       'RYYN'::dna AS dna3;

-- ======================================================
-- Equality tests
-- ======================================================
-- SELECT equals('ACGT'::dna, 'ACGT'::dna) AS equality_expected_true;
-- SELECT equals('ACGT'::dna, 'AGGT'::dna) AS equality_expected_false;

-- -- Operator form (=)
-- SELECT 'ACGT'::dna = 'ACGT'::dna AS equality_op_true;
-- SELECT 'ACGT'::dna = 'AGGT'::dna AS equality_op_false;

-- ======================================================
-- Length tests
-- ======================================================
SELECT length('A'::dna) AS length_1;
SELECT length('ACGT'::dna) AS length_4;
SELECT length('ACGTGGC'::dna) AS length_7;

-- ======================================================
-- Operator tests
-- ======================================================
-- The ^@ operator = “starts_with” if you defined it for DNA
-- SELECT 'ACGT'::dna ^@ 'A'::dna AS starts_with_op_true;
-- SELECT 'ACGT'::dna ^@ 'CG'::dna AS starts_with_op_false;

-- ======================================================
-- Table tests
-- ======================================================
DROP TABLE IF EXISTS dna_sequences;
CREATE TABLE dna_sequences (id serial, seq dna);

INSERT INTO dna_sequences (seq)
VALUES ('ACGT'::dna),
       ('AGGT'::dna),
       ('RYYN'::dna);

SELECT * FROM dna_sequences;

-- Query test
SELECT * FROM dna_sequences WHERE seq = 'ACGT'::dna;
