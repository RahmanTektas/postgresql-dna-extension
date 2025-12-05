DROP EXTENSION IF EXISTS dna_sequence CASCADE;
CREATE EXTENSION dna_sequence;
SET client_min_messages TO 'notice';

-- 1. Setup Environment
DROP TABLE IF EXISTS kmer_test;
CREATE TABLE kmer_test (val kmer);


-- 2. Insert Data


INSERT INTO kmer_test VALUES 
    ('AAAAA');

-- We insert specific patterns to test the Trie structure
-- Group 1: Poly-A variations (tests depth)
-- INSERT INTO kmer_test VALUES 
--     ('AAAAA'), ('AAAAT'), ('AAAAG'), ('AAAAC');

-- Group 2: Prefix splitting
-- These share 'TGC' but diverge at the 4th base
-- INSERT INTO kmer_test VALUES 
--     ('TGCA'), ('TGCT'), ('TGCG'), ('TGCC');

-- Group 3: Random noise to fill pages
-- INSERT INTO kmer_test VALUES 
--     ('CCGG'), ('GGCC'), ('ATAT'), ('TATA');

-- 3. Create the SP-GiST Index
CREATE INDEX idx_kmer_spgist
ON kmer_test
USING spgist (val kmer_spgist_ops);


-- CREATE INDEX idx_kmer_spgist
-- ON kmer_test
-- USING spgist (val public.kmer_spgist_ops);

-- 4. ANALYZE to update statistics
ANALYZE kmer_test;

-- 5. FORCE INDEX USAGE
-- This is crucial. Without this, PG might just scan the table because it's small.
SET enable_seqscan = OFF;

SELECT * FROM kmer_test ORDER BY val;  -- Should show all 12 rows
---------------------------------------------------------------------
-- TEST 1: Equality (=)
-- Should return exactly one row 'TGCA'
---------------------------------------------------------------------
SELECT 'Test 1: Equality' as test_name;
EXPLAIN (COSTS OFF) SELECT * FROM kmer_test WHERE val = 'TGCA'::kmer;
SELECT * FROM kmer_test WHERE val = 'TGCA'::kmer;

