EXTENSION   = kmer
MODULES 	= kmer
DATA        = kmer--1.0.sql kmer.control

PG_CONFIG ?= pg_config
PGXS = $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

reinstall:		# make reinstall to rebuild and install the extension
	make clean
	make
	sudo make install


test:		# make test to run the tests
	sudo -u postgres psql -d test_dna -c "DROP EXTENSION IF EXISTS kmer CASCADE;"
	sudo -u postgres psql -d test_dna -c "CREATE EXTENSION kmer;"
	sudo -u postgres psql -d test_dna -f tests/test_kmer.sql

run :
	make reinstall
	make test