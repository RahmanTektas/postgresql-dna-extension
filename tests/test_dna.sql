-- 1. CLEAN STATE (ONLY HERE)
DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;

\echo '========================================='
\echo '===          TEST: DNA TYPE           ==='
\echo '========================================='

-- ======================================================
-- PART 1: MANUAL UNIT TESTS
-- ======================================================

-- 1. Parsing & Output
\echo '\n--- Test 1: DNA parsing and output ---'
SELECT 'ACGT'::dna AS dna1,
       'ACGN'::dna AS dna2,
       'RYYN'::dna AS dna3;


-- 2. Case Sensitivity
\echo '\n--- Test 2: Case insensitivity (lower to upper) ---'
SELECT 'acgt'::dna AS lower_to_upper;


-- 3. Length Tests
\echo '\n--- Test 3: Length verification (Generated Column) ---'
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


\echo '\n--- Test 4: Additional length insertions ---'
INSERT INTO dna_with_length (seq)
VALUES
    ('ACGT'::dna),
    ('ACGTTT'::dna);

-- Display
SELECT * FROM dna_with_length;
TRUNCATE TABLE dna_with_length RESTART IDENTITY;


\echo '\n--- Test 5: Length with mixed case inputs ---'
INSERT INTO dna_with_length (seq)
VALUES
    ('acgta'::dna),
    ('acgtttta'::dna);

-- Display
SELECT * FROM dna_with_length;
TRUNCATE dna_with_length RESTART IDENTITY;


\echo '\n--- Test 6: Basic DNA table insertions ---'
DROP TABLE IF EXISTS dna_sequences;
CREATE TABLE dna_sequences (id serial, seq dna);

INSERT INTO dna_sequences (seq)
VALUES ('ACGT'::dna),
       ('AGGT'::dna),
       ('RYYN'::dna);

SELECT * FROM dna_sequences;


\echo '\n--- Test 7: DNA table with validation ---'
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


\echo '\n--- Test 8: DNA casting to text validation ---'
SELECT seq::text, seq::text = expected AS matches
FROM (
    SELECT seq,
           CASE WHEN id = 1 THEN 'ACGTACGT'
                WHEN id = 2 THEN 'TTGCA'
                WHEN id = 3 THEN 'GGGAAA'
           END AS expected
    FROM test_dna
) t;


-- ======================================================
-- PART 2: DATASET TESTS
-- ======================================================

\echo '\n--- Test 9: Loading dataset from file ---'
DROP TABLE IF EXISTS test_dna;

CREATE TABLE test_dna(
    id bigserial PRIMARY KEY,
    seq dna
);

-- Note: Ensure the file exists at this path in your docker container/environment
COPY test_dna(seq)
FROM '/extension/dna_sequence/scriptpy/kmer_36.txt'
WITH (FORMAT text);


\echo '\n--- Test 10: Select random sample from table ---'
SELECT seq
FROM test_dna
ORDER BY random()
LIMIT 40;


\echo '\n--- Test 11: Select length from random sample ---'
SELECT length(seq)
FROM test_dna
ORDER BY random()
LIMIT 1;


\echo '\n--- Test 12: Integration - Insert random samples into length table ---'
INSERT INTO dna_with_length (seq)
SELECT seq
FROM test_dna
ORDER BY random()
LIMIT 10;

SELECT * FROM dna_with_length;
