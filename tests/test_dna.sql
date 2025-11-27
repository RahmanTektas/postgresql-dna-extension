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


\echo 'test 2 lower to upper'
SELECT 'acgt'::dna AS lower_to_upper;

-- ======================================================
-- Length tests
-- ======================================================
\echo 'test 3 test lenght'

-- Drop correct table name
DROP TABLE IF EXISTS dna_with_length;

-- Create table with generated column
CREATE TABLE dna_with_length (
    id   bigserial PRIMARY KEY,
    seq  dna NOT NULL,
    len  integer GENERATED ALWAYS AS (length(seq)) STORED
);

-- Insert rows
INSERT INTO dna_with_length (seq)
VALUES
    ('ACGT'::dna),
    ('ACGTTT'::dna),
    ('NNNN'::dna);

-- Display
SELECT * FROM dna_with_length;


\echo 'test 4 test lenght'
INSERT INTO dna_with_length (seq)
VALUES
    ('ACGT'::dna),
    ('ACGTTT'::dna);

-- Display
SELECT * FROM dna_with_length;
TRUNCATE TABLE dna_with_length RESTART IDENTITY;


\echo 'test 5 test lenght with lower char'
INSERT INTO dna_with_length (seq)
VALUES
    ('acgta'::dna),
    ('acgtttta'::dna);

-- Display
SELECT * FROM dna_with_length;
TRUNCATE dna_with_length RESTART IDENTITY;


-- ======================================================
-- Operator tests
-- ======================================================
-- The ^@ operator = “starts_with” if you defined it for DNA
-- SELECT 'ACGT'::dna ^@ 'A'::dna AS starts_with_op_true;
-- SELECT 'ACGT'::dna ^@ 'CG'::dna AS starts_with_op_false;

-- ======================================================
-- Table tests
-- ======================================================
\echo '--- Test 6: DNA in table ---'
DROP TABLE IF EXISTS dna_sequences;
CREATE TABLE dna_sequences (id serial, seq dna);

INSERT INTO dna_sequences (seq)
VALUES ('ACGT'::dna),
       ('AGGT'::dna),
       ('RYYN'::dna);

SELECT * FROM dna_sequences;

\echo '--- Test 7: DNA table with validation ---'
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

\echo '--- Test 8: DNA casting to text validation ---'
SELECT seq::text, seq::text = expected AS matches
FROM (
    SELECT seq,
           CASE WHEN id = 1 THEN 'ACGTACGT'
                WHEN id = 2 THEN 'TTGCA'
                WHEN id = 3 THEN 'GGGAAA'
           END AS expected
    FROM test_dna
) t;


\echo '--- Test 9: test with dataset ---'

DROP TABLE IF EXISTS test_dna;

CREATE TABLE test_dna(
    id bigserial PRIMARY KEY,
    seq dna
);

COPY test_dna(seq)
FROM '/extension/dna_sequence/scriptpy/kmers_36.txt'
WITH (FORMAT text);

\echo '--- Test 10: select a random sample dna sequence from table ---'

SELECT seq
FROM test_dna
ORDER BY random()
LIMIT 40;


\echo '--- Test 10: select a random sample dna sequence from table ---'

SELECT length(seq)
FROM test_dna
ORDER BY random()
LIMIT 1;

\echo '--- Test 11: select size from dataset ---'


INSERT INTO dna_with_length (seq)
SELECT seq
FROM test_dna
ORDER BY random()
LIMIT 10;


SELECT * 
FROM dna_with_length;

TRUNCATE dna_with_length RESTART IDENTITY;


