#!/bin/bash

# Colors for better visibility
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Configuration
CONTAINER="info-h417-dna-extension"
USER="postgres"

# Function to wait for user input
wait_user() {
    echo -e "${YELLOW}>>> Press ENTER to continue...${NC}"
    read -r
}

# Function to execute an SQL command
run_sql() {
    echo -e "${GREEN}$1${NC}"
    # Using || true to prevent script exit on expected SQL errors
    docker exec -u $USER $CONTAINER psql -c "$2" || true
}

clear
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║              INTERACTIVE DEMO: DNA Sequence Extension              ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""
wait_user

# =========================================================================
# RESET
# =========================================================================
clear
echo -e "${BLUE}→ Initializing the extension...${NC}"
docker exec -u $USER $CONTAINER psql -c "DROP EXTENSION IF EXISTS dna_sequence CASCADE;"
docker exec -u $USER $CONTAINER psql -c "CREATE EXTENSION dna_sequence;"
echo ""
wait_user

# =========================================================================
# SECTION 1: DNA TYPE
# =========================================================================
clear
echo "┌─────────────────────────────────────────┐"
echo "│  SECTION 1: DNA TYPE                    │"
echo "└─────────────────────────────────────────┘"
echo ""
wait_user

echo -e "${GREEN}→ Test 1: Parsing valid 'ACGT', invalid 'ACGN' (N), invalid 'RYYN' (R/Y):${NC}"
run_sql "SELECT 'ACGT'::dna AS dna1, 'ACGN'::dna AS dna2, 'RYYN'::dna AS dna3;" "SELECT 'ACGT'::dna AS dna1, 'ACGN'::dna AS dna2, 'RYYN'::dna AS dna3;"
echo ""
wait_user

echo -e "${GREEN}→ Test 2: Case insensitivity (input 'acgt'):${NC}"
run_sql "SELECT 'acgt'::dna AS lower_to_upper;" "SELECT 'acgt'::dna AS lower_to_upper;"
echo ""
wait_user

echo -e "${GREEN}→ Test 3: Length of 'ACGTACGTACGT' (12 chars):${NC}"
run_sql "SELECT length('ACGTACGTACGT'::dna) AS length;" "SELECT length('ACGTACGTACGT'::dna) AS length;"
echo ""
wait_user

echo -e "${GREEN}→ Test 3-5: Table Storage, Generated Column & Mixed Case:${NC}"
run_sql "Creating table 'dna_with_length'..." "DROP TABLE IF EXISTS dna_with_length; CREATE TABLE dna_with_length (id bigserial PRIMARY KEY, seq dna NOT NULL, len integer GENERATED ALWAYS AS (length(seq)) STORED);"
run_sql "Inserting 'ACGT', 'ACGTTT', 'NNNN' (Error)..." "INSERT INTO dna_with_length (seq) VALUES ('ACGT'::dna), ('ACGTTT'::dna), ('NNNN'::dna);"
run_sql "Displaying table..." "SELECT * FROM dna_with_length;"
run_sql "Insertion mixed case..." "TRUNCATE TABLE dna_with_length RESTART IDENTITY; INSERT INTO dna_with_length (seq) VALUES ('acgta'::dna), ('acgtttta'::dna); SELECT * FROM dna_with_length;"
echo ""
wait_user

echo -e "${GREEN}→ Test 6: Basic DNA table insertions:${NC}"
run_sql "Inserting values..." "DROP TABLE IF EXISTS dna_sequences; CREATE TABLE dna_sequences (id serial, seq dna); INSERT INTO dna_sequences (seq) VALUES ('ACGT'::dna), ('AGGT'::dna), ('RYYN'::dna); SELECT * FROM dna_sequences;"
echo ""
wait_user

echo -e "${GREEN}→ Test 7: DNA table with validation:${NC}"
run_sql "Setup test_dna..." "DROP TABLE IF EXISTS test_dna; CREATE TABLE test_dna (id serial PRIMARY KEY, seq dna);"
run_sql "Inserting & Selecting..." "INSERT INTO test_dna (seq) VALUES ('ACGTACGT'::dna), ('TTGCA'::dna), ('GGGAAA'::dna); SELECT id, seq FROM test_dna;"
echo ""
wait_user

# =========================================================================
# SECTION 2: KMER TYPE
# =========================================================================
clear
echo "┌─────────────────────────────────────────┐"
echo "│  SECTION 2: KMER TYPE & OPERATORS       │"
echo "└─────────────────────────────────────────┘"
echo ""
wait_user

echo -e "${GREEN}→ Test 1: Parsing 'ACGT', 'aaaa', 'AGCT':${NC}"
run_sql "Parsing specific kmers..." "SELECT 'ACGT'::kmer AS k1, 'aaaa'::kmer AS k2, 'AGCT'::kmer AS k3;"
echo ""
wait_user

echo -e "${GREEN}→ Test 2: Equality function and operator (=):${NC}"
run_sql "Comparing using equals() and = ..." "SELECT equals('ACGT'::kmer, 'ACGT'::kmer) AS equality_expected_true, equals('ACGT'::kmer, 'AGCT'::kmer) AS equality_expected_false;"
run_sql "Operator check..." "SELECT 'ACGT'::kmer = 'ACGT'::kmer AS equality_operator_true, 'ACGT'::kmer = 'AGCT'::kmer AS equality_operator_false;"
echo ""
wait_user

echo -e "${GREEN}→ Test 3: Length verification:${NC}"
run_sql "Checking lengths..." "SELECT length('ACGT'::kmer) AS expected_length_of_kmer_4, length('A'::kmer) AS expected_length_of_kmer_1;"
echo ""
wait_user

echo -e "${GREEN}→ Test 4: Starts_with function and operator (^@):${NC}"
run_sql "Checking starts_with()..." "SELECT starts_with('ACGT'::kmer, 'A'::kmer) AS starts_with_expected_true, starts_with('ACGT'::kmer, 'ACGTA'::kmer) AS starts_with_expected_false;"
run_sql "Checking operator ^@ ..." "SELECT 'ACGT'::kmer ^@ 'AC'::kmer AS starts_with_op_test1_true, 'ACGT'::kmer ^@ 'CG'::kmer AS starts_with_op_test2_false;"
echo ""
wait_user

echo -e "${GREEN}→ Test 5: Basic Table Operations:${NC}"
run_sql "Inserting specific kmers..." "DROP TABLE IF EXISTS kmers_manual; CREATE TABLE kmers_manual (id serial, seq kmer); INSERT INTO kmers_manual (seq) VALUES ('ACGT'::kmer), ('AAAA'::kmer), ('AGCT'::kmer);"
run_sql "Selecting 'ACGT'..." "SELECT * FROM kmers_manual WHERE seq = 'ACGT'::kmer;"
echo ""
wait_user

echo -e "${GREEN}→ Test 6: generate_kmers(k=3) functionality:${NC}"
run_sql "Generating 3-mers from 'ACGTACGT'..." "SELECT k.kmer FROM generate_kmers('ACGTACGT'::dna, 3) AS k(kmer);"
echo ""
wait_user

echo -e "${GREEN}→ Test 7: Grouping and Counting (GROUP BY):${NC}"
run_sql "Counting occurrences in 'ACGTACGT'..." "WITH kmer_counts AS (SELECT k.kmer, count(*) as count FROM generate_kmers('ACGTACGT'::dna, 3) AS k(kmer) GROUP BY k.kmer) SELECT * FROM kmer_counts ORDER BY count DESC;"
echo ""
wait_user


# =========================================================================
# SECTION 3: QKMER TYPE
# =========================================================================
clear
echo "┌─────────────────────────────────────────┐"
echo "│  SECTION 3: QKMER TYPE (WILDCARDS)      │"
echo "└─────────────────────────────────────────┘"
echo ""
wait_user

echo -e "${GREEN}→ Test 1: Parsing 'ANGT', 'ARGT', 'AYGT' as qkmer:${NC}"
run_sql "Parsing specific qkmers..." "SELECT 'ANGT'::qkmer AS qk1, 'ARGT'::qkmer AS qk2, 'AYGT'::qkmer AS qk3;"
echo ""
wait_user

echo -e "${GREEN}→ Test 2: Length of qkmer 'AYGT':${NC}"
run_sql "Checking length..." "SELECT length('AYGT'::qkmer) AS expected_length_of_qkmer_4;"
echo ""
wait_user

echo -e "${GREEN}→ Test 3: IUPAC Codes: 'ANGT' (N), 'ARGT' (R), 'AYGT' (Y):${NC}"
run_sql "Parsing wildcards..." "SELECT 'ANGT'::qkmer AS wild_N, 'ARGT'::qkmer AS wild_R, 'AYGT'::qkmer AS wild_Y;"
echo ""
wait_user

echo -e "${GREEN}→ Test 4: Contains Operator (@>):${NC}"
run_sql "Checking coverage..." "SELECT 'ANGT'::qkmer @> 'ACGT'::kmer AS match_N_true, 'ARGT'::qkmer @> 'AGGT'::kmer AS match_R_true, 'ARGT'::qkmer @> 'AGGTT'::kmer AS match_length_false, 'ARGT'::qkmer @> 'ATGT'::kmer AS match_R_false;"
echo ""
wait_user

echo -e "${GREEN}→ Test 5: Multiple ambiguous positions:${NC}"
run_sql "Complex matching..." "SELECT 'RNYT'::qkmer @> 'AACT'::kmer AS true_test_same_len_1, 'RNYT'::qkmer @> 'GGTT'::kmer AS true_test_same_len_2, 'RNYT'::qkmer @> 'CCCT'::kmer AS false_test_wrong_first_base, 'AANN'::qkmer @> 'AAGT'::kmer AS true_test_N_any, 'AARY'::qkmer @> 'AAGT'::kmer AS true_test_RY, 'AARY'::qkmer @> 'AACT'::kmer AS false_test_RY_last_base;"
echo ""
wait_user

echo -e "${GREEN}→ Test 6: Full wildcard N-pattern:${NC}"
run_sql "Full N wildcard..." "SELECT 'NNNN'::qkmer @> 'ACGT'::kmer AS full_N_true, 'NNNN'::qkmer @> 'TTTT'::kmer AS full_N_true2;"
echo ""
wait_user


# =========================================================================
# SECTION 4: DATASET INTEGRATION
# =========================================================================
clear
echo "┌─────────────────────────────────────────┐"
echo "│  SECTION 4: DATASET (DNA, KMER, QKMER)  │"
echo "└─────────────────────────────────────────┘"
echo ""
wait_user

# --- PART 1: DNA DATASET ---
echo -e "${BLUE}=== PART 1: DNA DATASET ===${NC}"

echo -e "${GREEN}→ Test 9: Loading dataset from file:${NC}"
run_sql "Creating table 'test_dna' & COPY from 'dna_36.txt'..." "DROP TABLE IF EXISTS test_dna; CREATE TABLE test_dna(id bigserial PRIMARY KEY, seq dna); COPY test_dna(seq) FROM '/extension/scriptpy/dna_36.txt' WITH (FORMAT text);"
echo ""
wait_user

echo -e "${GREEN}→ Test 10: Select random sample from table:${NC}"
run_sql "Selecting 5 random sequences..." "SELECT seq FROM test_dna ORDER BY random() LIMIT 5;"
echo ""
wait_user

echo -e "${GREEN}→ Test 11: Select length from random sample:${NC}"
run_sql "Checking length of 1 random sequence..." "SELECT length(seq) FROM test_dna ORDER BY random() LIMIT 1;"
echo ""
wait_user

echo -e "${GREEN}→ Test 12: Integration - Insert random samples into length table:${NC}"
run_sql "Setup: Ensuring dna_with_length exists..." "TRUNCATE TABLE dna_with_length RESTART IDENTITY;"
run_sql "Inserting & verifying..." "INSERT INTO dna_with_length (seq) SELECT seq FROM test_dna ORDER BY random() LIMIT 5; SELECT * FROM dna_with_length;"
echo ""
wait_user

echo -e "${GREEN}→ Test 13: Empty sequence '' (Expecting Error):${NC}"
echo -e "${RED}Note: DNA constraints forbid empty sequences.${NC}"
run_sql "Attempting to cast empty string to DNA..." "SELECT ''::dna;"
echo ""
wait_user

# --- PART 2: KMER GENERATION ---
clear
echo -e "${BLUE}=== PART 2: KMER GENERATION LOGIC ===${NC}"

echo -e "${GREEN}→ Test 8, 9, 10: Manual Generation with k=3, k=4, k=32:${NC}"
run_sql "Setup..." "DROP TABLE IF EXISTS kmer_dataset; CREATE TABLE kmer_dataset (id serial PRIMARY KEY, val kmer);"
run_sql "Inserting k=3 ('AAAAAA' -> 'AAA')..." "INSERT INTO kmer_dataset (val) SELECT k.kmer FROM generate_kmers('AAAAAA'::dna, 3) AS k(kmer);"
run_sql "Inserting k=4 ('ACGTACGT')..." "INSERT INTO kmer_dataset (val) SELECT k.kmer FROM generate_kmers('ACGTACGT'::dna, 4) AS k(kmer);"
run_sql "Inserting k=32 (max length)..." "INSERT INTO kmer_dataset (val) SELECT k.kmer FROM generate_kmers('ACGTACGTACGTACGTACGTACGTACGTACGT'::dna, 32) AS k(kmer);"
echo ""

echo -e "${GREEN}→ Test 11: Dataset Statistics (Expect 10 rows):${NC}"
run_sql "Checking min/max/count..." "SELECT count(*) as total_rows_should_be_10, min(length(val)) as min_len_should_be_3, max(length(val)) as max_len_should_be_32 FROM kmer_dataset;"
echo ""
wait_user

echo -e "${GREEN}→ Test 12: Filtering generated data starts_with 'AAA':${NC}"
run_sql "Counting 'AAA'..." "SELECT count(*) as count_starts_with_AAA_should_be_4 FROM kmer_dataset WHERE val ^@ 'AAA'::kmer;"
echo ""
wait_user

echo -e "${GREEN}→ Test 13: Grouping generated duplicates 'AAA':${NC}"
run_sql "Grouping..." "SELECT val, count(*) FROM kmer_dataset WHERE val = 'AAA'::kmer GROUP BY val;"
echo ""
wait_user

echo -e "${BLUE}=== PART 3: KMER DATASET FROM FILE ===${NC}"

echo -e "${GREEN}→ Test 14 (Setup): Reloading DNA Source File:${NC}"
run_sql "Resetting test_dna..." "DROP TABLE IF EXISTS test_dna; CREATE TABLE test_dna(id bigserial PRIMARY KEY, seq dna); COPY test_dna(seq) FROM '/extension/scriptpy/dna_36.txt' WITH (FORMAT text);"
echo ""

echo -e "${GREEN}→ Test 15: Creating empty k-mer dataset table:${NC}"
run_sql "Creating other_kmer_dataset..." "DROP TABLE IF EXISTS other_kmer_dataset; CREATE TABLE other_kmer_dataset(id bigserial PRIMARY KEY, seq kmer);"
echo ""
wait_user

echo -e "${GREEN}→ Test 16 & 17: Generating and inspecting 4-mers (Limit 1000):${NC}"
run_sql "Inserting 4-mers..." "INSERT INTO other_kmer_dataset (seq) SELECT k.kmer FROM test_dna d, generate_kmers(d.seq, 4) AS k(kmer) LIMIT 1000;"
run_sql "Verifying 4-mers..." "SELECT * FROM other_kmer_dataset WHERE length(seq) = 4 LIMIT 5;"
echo ""
wait_user

echo -e "${GREEN}→ Test 18 & 19: Generating and inspecting 10-mers (Limit 1000):${NC}"
run_sql "Inserting 10-mers..." "INSERT INTO other_kmer_dataset (seq) SELECT k.kmer FROM test_dna d, generate_kmers(d.seq, 10) AS k(kmer) LIMIT 1000;"
run_sql "Verifying 10-mers..." "SELECT * FROM other_kmer_dataset WHERE length(seq) = 10 LIMIT 5;"
echo ""
wait_user

echo -e "${GREEN}→ Test 20 & 21: Generating and inspecting 19-mers (Limit 1000):${NC}"
run_sql "Inserting 19-mers..." "INSERT INTO other_kmer_dataset (seq) SELECT k.kmer FROM test_dna d, generate_kmers(d.seq, 19) AS k(kmer) LIMIT 1000;"
run_sql "Verifying 19-mers..." "SELECT * FROM other_kmer_dataset WHERE length(seq) = 19 LIMIT 5;"
echo ""
wait_user

echo -e "${GREEN}→ Test 22 & 23: Generating 33-mers (Should be empty):${NC}"
run_sql "Attempting to insert 33-mers..." "INSERT INTO other_kmer_dataset (seq) SELECT k.kmer FROM test_dna d, generate_kmers(d.seq, 33) AS k(kmer) LIMIT 1000;"
run_sql "Verifying table has no 33-mers..." "SELECT * FROM other_kmer_dataset WHERE length(seq) = 33 LIMIT 5;"
echo ""
wait_user

echo -e "${GREEN}→ Test 24: Counting all generated kmers grouped by length:${NC}"
run_sql "Distribution of Kmers..." "SELECT length(seq) AS kmer_length, count(*) AS count_should_be_1000 FROM other_kmer_dataset GROUP BY length(seq) ORDER BY kmer_length LIMIT 5;"
echo ""
wait_user

# --- PART 4: QKMER DATASET ---
clear
echo -e "${BLUE}=== PART 4: QKMER DATASET OPERATIONS ===${NC}"

echo -e "${GREEN}→ Test 7: Loading qkmer targets from file:${NC}"
run_sql "Loading qkmer_targets_dataset..." "DROP TABLE IF EXISTS qkmer_targets_dataset; CREATE TABLE qkmer_targets_dataset (id bigserial PRIMARY KEY, val qkmer); COPY qkmer_targets_dataset(val) FROM '/extension/scriptpy/qkmeroutput.txt' WITH (FORMAT text);"
echo ""
wait_user

echo -e "${GREEN}→ Test 8: Verify Dataset Loading:${NC}"
run_sql "Checking count..." "SELECT count(*) as total_rows_loaded FROM qkmer_targets_dataset;"
echo ""
wait_user

echo -e "${GREEN}→ Test 9: Search dataset using invalid 'AAAN' (Expecting Error):${NC}"
echo -e "${RED}Note: 'AAAN' is cast to kmer (exact), but N is not valid in Kmer type.${NC}"
run_sql "Executing invalid query..." "SELECT val AS qkmer_including_AAAN FROM qkmer_targets_dataset WHERE val @> 'AAAN'::kmer LIMIT 5;"
echo ""
wait_user

echo -e "${GREEN}→ Test 10: Check if 'CGTANY' covers 'CGTAAC':${NC}"
run_sql "Creating temp table with 'CGTANY' & checking..." "DROP TABLE IF EXISTS qkmer_with_cgtany; CREATE TABLE qkmer_with_cgtany (id bigserial PRIMARY KEY, val qkmer); INSERT INTO qkmer_with_cgtany (val) VALUES ('CGTANY'::qkmer); SELECT val AS qkmer_covering_CGTA FROM qkmer_with_cgtany WHERE val @> 'CGTAAC'::kmer;"
echo ""
wait_user

echo -e "${GREEN}→ Test 11: Dataset: Search qkmers covering 'CGTA':${NC}"
run_sql "Searching in dataset..." "SELECT val AS qkmer_covering_CGTA FROM qkmer_targets_dataset WHERE val @> 'CGTA'::kmer LIMIT 5;"
echo ""
wait_user

echo -e "${GREEN}→ Test 12: Dataset: Search qkmers covering 'TCGA':${NC}"
run_sql "Searching in dataset..." "SELECT val AS qkmer_covering_TCGA FROM qkmer_targets_dataset WHERE val @> 'TCGA'::kmer LIMIT 5;"
echo ""
wait_user


# =========================================================================
# SECTION 5: SP-GIST INDEX
# =========================================================================
clear
echo "┌─────────────────────────────────────────┐"
echo "│  SECTION 5: SP-GIST INDEX               │"
echo "└─────────────────────────────────────────┘"
echo ""
wait_user

echo -e "${GREEN}→ Setup: Disable Sequential Scan${NC}"
run_sql "SET enable_seqscan = OFF;" "SET enable_seqscan = OFF;"
echo ""

echo -e "${GREEN}→ Test 1: Setup and Operator Class Check:${NC}"
run_sql "Checking kmer_spgist_ops..." "DROP TABLE IF EXISTS kmer_test; CREATE TABLE kmer_test (val kmer); SELECT opcname AS operator_class_name, amname AS access_method, opcintype::regtype AS input_type, n.nspname AS schema FROM pg_opclass c JOIN pg_am a ON c.opcmethod = a.oid JOIN pg_namespace n ON c.opcnamespace = n.oid WHERE opcname = 'kmer_spgist_ops';"
echo ""
wait_user

echo -e "${GREEN}→ Test 2 & 3: Inserting Patterns & Creating Index:${NC}"
run_sql "Inserting patterns..." "INSERT INTO kmer_test VALUES ('AAAAA'), ('AAAAT'), ('AAAAG'), ('AAAAC'), ('TGCA'), ('TGCT'), ('TGCG'), ('TGCC'), ('CCGG'), ('GGCC'), ('ATAT'), ('TATA');"
run_sql "Creating Index..." "CREATE INDEX idx_kmer_spgist ON kmer_test USING spgist (val public.kmer_spgist_ops);"
echo ""
wait_user

echo -e "${GREEN}→ Test 4: Equality Search (=) 'TGCA':${NC}"
run_sql "Explain & Select..." "EXPLAIN (COSTS OFF) SELECT val AS matched_kmer FROM kmer_test WHERE val = 'TGCA'; SELECT val AS matched_kmer FROM kmer_test WHERE val = 'TGCA';"
echo ""
wait_user

echo -e "${GREEN}→ Test 5: Prefix Search (^@) 'TGC':${NC}"
run_sql "Explain & Select..." "EXPLAIN (COSTS OFF) SELECT val AS prefix_matched_kmer FROM kmer_test WHERE val ^@ 'TGC'; SELECT val AS prefix_matched_kmer FROM kmer_test WHERE val ^@ 'TGC';"
echo ""
wait_user

echo -e "${GREEN}→ Test 6 & 7: QKmer Search ('AAAAN', 'TGCR'):${NC}"
run_sql "Wildcard N..." "EXPLAIN (COSTS OFF) SELECT val AS wildcard_matched_kmer FROM kmer_test WHERE 'AAAAN'::qkmer @> val; SELECT val AS wildcard_matched_kmer FROM kmer_test WHERE 'AAAAN'::qkmer @> val;"
run_sql "Wildcard R..." "SELECT val AS wildcard_R_matched_kmer FROM kmer_test WHERE 'TGCR'::qkmer @> val;"
echo ""
wait_user

echo -e "${GREEN}→ Test 8: Length Mismatch (Expect 0 rows):${NC}"
run_sql "Searching 'AAA'::qkmer..." "SELECT val AS length_mismatch_kmer FROM kmer_test WHERE 'AAA'::qkmer @> val;"
echo ""
wait_user


echo -e "${BLUE}=== PART 2: SP-GIST DATASET TESTS ===${NC}"

echo -e "${GREEN}→ Test 9-11: Prepare Dataset (Load & Generate k=5):${NC}"
run_sql "Creating tables & Loading..." "DROP TABLE IF EXISTS spg_kmer_dataset; CREATE TABLE spg_kmer_dataset(id bigserial PRIMARY KEY, seq kmer); DROP TABLE IF EXISTS test_dna; CREATE TABLE test_dna(id bigserial PRIMARY KEY, seq dna); COPY test_dna(seq) FROM '/extension/scriptpy/dna_36.txt' WITH (FORMAT text);"
run_sql "Generating Kmers..." "INSERT INTO spg_kmer_dataset(seq) SELECT k.kmer AS kmer_value FROM test_dna d, generate_kmers(d.seq, 5) AS k(kmer);"
echo ""
wait_user

echo -e "${GREEN}→ Test 12: Create SP-GiST index for dataset:${NC}"
run_sql "Building Index..." "DROP INDEX IF EXISTS spg_kmer_dataset_idx; CREATE INDEX spg_kmer_dataset_idx ON spg_kmer_dataset USING spgist (seq);"
echo ""
wait_user

echo -e "${GREEN}→ Test 13: Prefix query 'ACG':${NC}"
run_sql "Explain Analyze..." "EXPLAIN ANALYZE SELECT id AS record_id, seq AS prefix_matched_kmer FROM spg_kmer_dataset WHERE seq ^@ 'ACG'::kmer LIMIT 10;"
run_sql "Result..." "SELECT DISTINCT seq AS prefix_matched_kmer FROM spg_kmer_dataset WHERE seq ^@ 'ACG'::kmer LIMIT 10;"
echo ""
wait_user

echo -e "${GREEN}→ Test 14: Random equality lookup:${NC}"
run_sql "Looking up random kmer..." "EXPLAIN ANALYZE SELECT id AS record_id, seq AS random_matched_kmer FROM spg_kmer_dataset WHERE seq = (SELECT seq FROM spg_kmer_dataset ORDER BY random() LIMIT 1);"
echo ""
wait_user

echo -e "${GREEN}→ Test 15: Qkmer pattern search 'RYNNN':${NC}"
run_sql "Explain Analyze..." "EXPLAIN ANALYZE SELECT id AS record_id, seq AS qkmer_pattern_match FROM spg_kmer_dataset WHERE seq <@ 'RYNNN'::qkmer LIMIT 10;"
run_sql "Result..." "SELECT DISTINCT seq AS qkmer_pattern_match FROM spg_kmer_dataset WHERE seq <@ 'RYNNN'::qkmer LIMIT 10;"
echo ""
wait_user

echo -e "${GREEN}→ Test 16: Count kmers grouped by length:${NC}"
run_sql "Grouping..." "SELECT length(seq) AS kmer_length, count(*) AS count_by_length FROM spg_kmer_dataset GROUP BY length(seq) ORDER BY kmer_length;"
echo ""
wait_user

echo -e "${GREEN}→ Test 17: Random sample of dataset:${NC}"
run_sql "Sampling..." "SELECT id AS record_id, seq AS random_sample_kmer FROM spg_kmer_dataset ORDER BY random() LIMIT 5;"
echo ""
wait_user


# =========================================================================
# END
# =========================================================================
clear
echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║                  DEMONSTRATION COMPLETED                           ║"
echo "╠════════════════════════════════════════════════════════════════════╣"
echo "║  ✓ DNA Type   : Validation and storage                             ║"
echo "║  ✓ KMER Type  : Generation and analysis                            ║"
echo "║  ✓ QKMER Type : Pattern matching with wildcards                    ║"
echo "║  ✓ DATASET    : Integrated tests (DNA, Kmer Gen, Qkmer Targets)    ║"
echo "║  ✓ Index      : SP-GIST for optimized searches                     ║"
echo "╚════════════════════════════════════════════════════════════════════╝"
echo ""