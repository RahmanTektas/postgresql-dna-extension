EXTENSION   = dna_sequence
MODULES 	= dna_sequence
DATA        = dna_sequence--1.0.sql dna_sequence.control

PG_CONFIG ?= pg_config
PGXS = $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

TEST_FOLDER=../tests/


reinstall:		# make reinstall to rebuild and install the extension
	make clean
	make
	sudo make install


test:		# make test to run the tests
	sudo -u postgres psql -d test_dna -c "DROP EXTENSION IF EXISTS dna_sequence CASCADE;"
	sudo -u postgres psql -d test_dna -c "CREATE EXTENSION dna_sequence;"
	sudo -u postgres psql -d test_dna -f $(TEST_FOLDER)test_kmer.sql
	sudo -u postgres psql -d test_dna -f $(TEST_FOLDER)test_dna.sql


run :
	make reinstall
	make test