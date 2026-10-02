# PostgreSQL DNA Extension

A native PostgreSQL extension written in C for storing and querying genomic sequences.

This was developed as a four-person ULB Database Systems Architecture project. The repository history preserves the individual contributions of the team members.

## Features

- custom variable-length `dna`, `kmer`, and `qkmer` PostgreSQL types
- input validation and C-level input/output functions
- sequence length, equality, prefix, containment, and k-mer generation operations
- SQL operators backed by native C functions
- hash and B-tree support for exact values
- SP-GiST indexing for prefix and pattern-oriented searches
- FASTQ preprocessing utilities
- SQL tests using synthetic and real genomic data

## Structure

```text
dna_sequence/
  dna_sequence.c          type representation and sequence operations
  spgist.c                SP-GiST callbacks
  dna_sequence--1.0.sql   SQL types, functions, operators and operator classes
  dna_sequence.h          shared C declarations

tests/                    SQL test suites
scriptpy/                 dataset preparation utilities
Makefile                  build, Docker and test commands
```

## Build

The original project used a course-provided Docker environment. In a compatible PostgreSQL development environment with PGXS and a C toolchain:

```bash
cd dna_sequence
make
sudo make install
psql -d YOUR_DATABASE -c "CREATE EXTENSION dna_sequence;"
```

For the Docker workflow configured in the repository:

```bash
make docker-reinstall
```

## Example

```sql
CREATE EXTENSION dna_sequence;

SELECT length('ACGTACGT'::dna);

SELECT *
FROM generate_kmers('ACGTACGT'::dna, 3);
```

## Tests

Run the complete suite:

```bash
make run
```

Focused targets are also available:

```bash
make test-dna
make test-kmer
make test-qkmer
make test-spgist
```

The project works directly with PostgreSQL extension APIs, varlena data representation, SQL operator classes, and index callbacks rather than implementing the functionality only at application level.
