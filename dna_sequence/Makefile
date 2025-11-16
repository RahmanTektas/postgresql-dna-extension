EXTENSION   = dna_sequence 
MODULES     = dna_sequence
DATA        = dna_sequence--1.0.sql dna_sequence.control

PG_CONFIG ?= pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs 2>/dev/null)
ifneq ($(PGXS),)
include $(PGXS)
endif

TEST_FOLDER    = ../tests

CONTAINER_NAME = info-h417-dna-extension

TEST_DB        = postgres

TEST_DOCKER_PSQL    = docker exec -i -u postgres $(CONTAINER_NAME) psql 

docker-reinstall:
	docker exec -u root $(CONTAINER_NAME) bash -lc "cd /extension/dna_sequence && make && make install"

test:
	$(TEST_DOCKER_PSQL) < $(TEST_FOLDER)/test_dna.sql
	$(TEST_DOCKER_PSQL) < $(TEST_FOLDER)/test_kmer.sql
	$(TEST_DOCKER_PSQL) < $(TEST_FOLDER)/test_qkmer.sql

# Rebuild + tests
run: docker-reinstall test
