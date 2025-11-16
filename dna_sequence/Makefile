EXTENSION   = dna_sequence
MODULES     = dna_sequence
DATA        = dna_sequence--1.0.sql dna_sequence.control

PG_CONFIG ?= pg_config
PGXS = $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

TEST_DB = postgres
TEST_FOLDER = ../tests/

# -----------------------------
# Build & Install (no sudo)
# -----------------------------
reinstall:
	make clean
	make
	make install

# -----------------------------
# Run tests inside Docker
# -----------------------------
test:
	psql -U postgres -d $(TEST_DB) -c "DROP EXTENSION IF EXISTS dna_sequence CASCADE;"
	psql -U postgres -d $(TEST_DB) -c "CREATE EXTENSION dna_sequence;"
	psql -U postgres -d $(TEST_DB) -f $(TEST_FOLDER)test_dna.sql
	psql -U postgres -d $(TEST_DB) -f $(TEST_FOLDER)test_kmer.sql
	psql -U postgres -d $(TEST_DB) -f $(TEST_FOLDER)test_qkmer.sql

# -----------------------------
# Shortcut: rebuild + test
# -----------------------------
run: reinstall test
