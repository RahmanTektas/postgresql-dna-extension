# PostgreSQL DNA Extension

A native PostgreSQL extension written in C for storing, transforming, and querying genomic sequences with domain-specific data types and index support.

> Academic team project for Database Systems Architecture at ULB. Developed by a four-person team; the repository history and contributor graph preserve individual authorship.

## Problem

Genomic sequences are often stored as plain text. That is convenient, but it prevents the database from understanding biological constraints or choosing specialized access paths. This extension moves that logic into PostgreSQL through native types, operators, functions, and an SP-GiST operator class.

## What the extension provides

- dna, kmer, and qkmer variable-length PostgreSQL types
- validation and C-level input/output functions
- sequence length, equality, prefix, containment, and k-mer generation operations
- SQL operators backed by native C functions
- hash and B-tree integration for exact values
- SP-GiST indexing for prefix and pattern-oriented searches
- utilities for transforming FASTQ data into loadable DNA and k-mer datasets
- automated SQL tests on synthetic and real genomic inputs

## Architecture

| Component | Responsibility |
| --- | --- |
| dna_sequence/dna_sequence.c | Type representation, parsing, validation, comparison, and sequence operations |
| dna_sequence/spgist.c | SP-GiST callbacks for indexed search |
| dna_sequence/dna_sequence--1.0.sql | SQL types, functions, operators, and operator classes |
| dna_sequence/dna_sequence.h | Shared C structures and declarations |
| tests/ | DNA, k-mer, qkmer, and SP-GiST test suites |
| scriptpy/ | FASTQ extraction and dataset-generation utilities |
| Makefile | Build, installation, Docker, and test automation |

## Engineering highlights

- PostgreSQL-compatible variable-length values using varlena layouts
- immutable, strict, and parallel-safe SQL functions where applicable
- domain validation at the type boundary rather than in application code
- IUPAC-aware query patterns through qkmer values
- index-aware operators integrated with PostgreSQL query execution
- reproducible command-line workflows for compilation and testing

## Build and install

The original project was developed in a course-provided Docker image. In a compatible PostgreSQL development environment with PGXS and a C toolchain:

    cd dna_sequence
    make
    sudo make install
    psql -d YOUR_DATABASE -c "CREATE EXTENSION dna_sequence;"

For the original Docker workflow configured by the repository:

    make docker-reinstall

## Example usage

    CREATE EXTENSION dna_sequence;

    SELECT length('ACGTACGT'::dna);

    SELECT *
    FROM generate_kmers('ACGTACGT'::dna, 3);

The complete SQL surface is declared in dna_sequence/dna_sequence--1.0.sql.

## Tests

Run the complete test suite:

    make run

Or run focused suites:

    make test-dna
    make test-kmer
    make test-qkmer
    make test-spgist

The test utilities support both controlled synthetic sequences and public sequencing data from the NCBI Sequence Read Archive.

## What this project demonstrates

This project goes beyond application-level SQL: it works with PostgreSQL's extension API, memory representation, operator semantics, planner-visible index support, C integration, and reproducible database testing.
