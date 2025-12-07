# Makefile for PostgreSQL DNA Sequence Extension

EXTENSION = dna_sequence
DATA = dna_sequence--1.0.sql
MODULE_big = dna_sequence
OBJS = dna_sequence.o spgist.o

PG_CONFIG = pg_config
PGXS := $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

# Docker Configuration
DOCKER_CONTAINER = info-h417-dna-extension
DOCKER_USER_ROOT = root
DOCKER_USER_POSTGRES = postgres
EXTENSION_DIR = /extension/dna_sequence
TESTS_DIR = ../tests
SCRIPTPY_DOCKER = /extension/

# Command to compile and install in Docker
.PHONY: docker-build
docker-build:
	docker exec -u $(DOCKER_USER_ROOT) $(DOCKER_CONTAINER) bash -lc "cd $(EXTENSION_DIR) && make clean && make && make install"

# Command to create the extension in PostgreSQL
.PHONY: docker-create-extension
docker-create-extension:
	docker exec -u $(DOCKER_USER_POSTGRES) $(DOCKER_CONTAINER) psql -c "DROP EXTENSION IF EXISTS $(EXTENSION) CASCADE; CREATE EXTENSION $(EXTENSION);"

# Command to compile, install AND create the extension (complete solution)
.PHONY: docker-reinstall
docker-reinstall: docker-build docker-sync-scriptpy docker-create-extension
	@echo "Extension compiled, installed, and created successfully!"

# Command to run all tests
.PHONY: run
run: docker-reinstall docker-sync-scriptpy
	@echo "Running DNA tests..."
	docker exec -i -u $(DOCKER_USER_POSTGRES) $(DOCKER_CONTAINER) psql < $(TESTS_DIR)/test_dna.sql
	@echo "\nRunning KMER tests..."
	docker exec -i -u $(DOCKER_USER_POSTGRES) $(DOCKER_CONTAINER) psql < $(TESTS_DIR)/test_kmer.sql
	@echo "\nRunning QKMER tests..."
	docker exec -i -u $(DOCKER_USER_POSTGRES) $(DOCKER_CONTAINER) psql < $(TESTS_DIR)/test_qkmer.sql
	@echo "\nRunning SPGIST tests..."
	docker exec -i -u $(DOCKER_USER_POSTGRES) $(DOCKER_CONTAINER) psql < $(TESTS_DIR)/test_spgist.sql


.PHONY: docker-sync-scriptpy
docker-sync-scriptpy:
	@echo "Syncing scriptpy/ into Docker container..."
	docker exec -u $(DOCKER_USER_ROOT) $(DOCKER_CONTAINER) mkdir -p $(SCRIPTPY_DOCKER)
	docker cp ../scriptpy/ $(DOCKER_CONTAINER):$(SCRIPTPY_DOCKER)/
	@echo "✔ scriptpy/ synced to $(SCRIPTPY_DOCKER)"


# Command to test DNA only
.PHONY: test-dna
test-dna:
	docker exec -i -u $(DOCKER_USER_POSTGRES) $(DOCKER_CONTAINER) psql < $(TESTS_DIR)/test_dna.sql

# Command to test KMER only
.PHONY: test-kmer
test-kmer:
	docker exec -i -u $(DOCKER_USER_POSTGRES) $(DOCKER_CONTAINER) psql < $(TESTS_DIR)/test_kmer.sql

# Command to test QKMER only
.PHONY: test-qkmer
test-qkmer:
	docker exec -i -u $(DOCKER_USER_POSTGRES) $(DOCKER_CONTAINER) psql < $(TESTS_DIR)/test_qkmer.sql

# Command to test spgist only
.PHONY: test-spgist
test-spgist:
	docker exec -i -u $(DOCKER_USER_POSTGRES) $(DOCKER_CONTAINER) psql < $(TESTS_DIR)/test_spgist.sql

# Command to connect to psql
.PHONY: psql
psql:
	docker exec -it $(DOCKER_CONTAINER) psql -U $(DOCKER_USER_POSTGRES)

# Command to clean everything
.PHONY: docker-clean
docker-clean:
	docker exec -u $(DOCKER_USER_ROOT) $(DOCKER_CONTAINER) bash -lc "cd $(EXTENSION_DIR) && make clean && rm -f *.o *.so *.bc"

# Help
.PHONY: help
help:
	@echo "Available commands:"
	@echo "  make docker-reinstall  - Compile, install AND create extension in PostgreSQL"
	@echo "  make run               - Reinstall extension and run all tests"
	@echo "  make test-dna          - Run DNA tests only"
	@echo "  make test-kmer         - Run KMER tests only"
	@echo "  make test-qkmer        - Run QKMER tests only"
	@echo "  make test-spgist       - Run spgist tests only"
	@echo "  make psql              - Connect to PostgreSQL"
	@echo "  make docker-clean      - Clean all compiled files"
	@echo "  make help              - Display this help message"
