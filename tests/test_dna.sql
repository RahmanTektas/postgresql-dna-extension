-- Clean state
DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;


\echo '========================================='
\echo '=== DNA TESTS ==========================='
\echo '========================================='

-- Basic parsing and output tests
\echo '--- Test 1: DNA parsing and output ---'
SELECT 'ACGT'::dna AS dna1,
       'ACGN'::dna AS dna2,
       'RYYN'::dna AS dna3;


-- ======================================================
-- Length tests
-- ======================================================
\echo '--- Test 2: DNA length function ---'
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
\echo '--- Test 3: DNA in table ---'
DROP TABLE IF EXISTS dna_sequences;
CREATE TABLE dna_sequences (id serial, seq dna);

INSERT INTO dna_sequences (seq)
VALUES ('ACGT'::dna),
       ('AGGT'::dna),
       ('RYYN'::dna);

SELECT * FROM dna_sequences;

\echo '--- Test 4: DNA table with validation ---'
DROP TABLE IF EXISTS test_dna;

CREATE TABLE test_dna (
    id serial PRIMARY KEY,
    seq dna
);

INSERT INTO test_dna (seq) VALUES
('ACGTACGT'::dna),
('TTGCA'::dna),
('GGGAAA'::dna);

-- Select to see if output is correct
SELECT id, seq FROM test_dna;

\echo '--- Test 5: DNA casting to text validation ---'
SELECT seq::text, seq::text = expected AS matches
FROM (
    SELECT seq,
           CASE WHEN id = 1 THEN 'ACGTACGT'
                WHEN id = 2 THEN 'TTGCA'
                WHEN id = 3 THEN 'GGGAAA'
           END AS expected
    FROM test_dna
) t;