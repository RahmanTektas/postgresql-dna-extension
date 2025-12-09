# postgresql-dna-extension
PostgreSQL extension implementing custom DNA types (dna, kmer, qkmer) for efficient genomic sequence representation and manipulation.

## Makefile Usage Guide
This extension includes a Makefile designed to simplify development inside the Docker container used for INFO-F417.
It provides commands to compile, install, reload, test, and clean the PostgreSQL DNA extension directly inside Docker.

All commands assume that the Docker container name is:
dna_container
and that the extension code is mounted into:
/extension/dna_sequence

## 1. Reinstall the Extension (Compile + Install + Create)
Recompiles the extension, installs it into PostgreSQL, syncs the Python scripts, and recreates the extension:

- ```make docker-reinstall```

## Equivalent to manually doing:

- ```make clean```
- ```make```
- ```make install```

- ```DROP EXTENSION …; CREATE EXTENSION …;```

## 2. Run All Tests

Runs DNA, KMER, QKMER, and SP-GiST tests in order:
- ```make run```

This is the recommended command to validate all parts of your implementation.


## 3. Run Individual Test Suites
DNA tests only:

- ```make test-dna```

KMER tests only:

- ```make test-kmer```

QKMER tests only:

- ```make test-qkmer```

SP-GiST index tests only:

- ```make test-spgist```


## 4. Sync Python Scripts into Docker
Copies the entire scriptpy/ folder into the Docker container:

- ```make docker-sync-scriptpy```

Useful when modifying:

- ```FASTQ → DNA extraction scripts```
- ```FASTQ → k-mer / q-mer generation tools```


## 5. Connect to PostgreSQL Inside Docker
- ```make psql```

This opens an interactive psql session as the postgres user inside the container.


## 6. Clean Compiled Files
- ```make docker-clean```

Removes all compiled artifacts (.o, .so, .bc) from the Docker extension directory.


## 7. Full Command List
- ```make help```

Displays a summary of all available commands.

# Build
## To do once
- install container info-h417-course-image:latest

- path_to_dna_sequence
- container_name

```bash
docker run -it \
  --name container_name \
  -v $(path_to_dna_sequence):/extension/dna_sequence \
info-h417-course-image:latest
```

## To run the extension

- ```docker start container_name```
- ```docker exec -it container_name service postgresql start```
- ```make run```


## 

# Testing with Synthetic and Real Genomic Data
To properly validate the PostgreSQL DNA extension, you should test it using both synthetic data and real-world sequencing data. This ensures that the operators, types, and SP-GiST index perform correctly across a wide range of scenarios and scales.

## 1. Synthetic DNA Sequences

For development and early debugging, you can easily generate random DNA sequences of arbitrary length. Synthetic data is useful for:
* verifying basic functionality of kmer and qkmer
* testing boundary cases (length, invalid bases, wildcard coverage)
* validating SP-GiST navigation on controlled patterns

## 2. Real Genomic Data from NCBI SRA
To ensure robustness on real datasets, test the extension using public sequencing reads from the NCBI Sequence Read Archive (SRA). These datasets reflect the noise, scale, and biological complexity encountered in actual genomic workloads.

### Install SRA Toolkit (Prerequisite)
If running inside the Docker container, run these commands first:

- ```apt-get update```
- ```apt-get install sra-toolkit```

### Download SRA data
First, create a specific data directory to keep the project clean:

- ```mkdir -p data```

Then, prefetch the data into this directory:

- ```prefetch SRR026760 --output-directory data/```

### Convert .sra → .fastq
- ```fasterq-dump data/SRR026760/SRR026760.sra -O data/```

### Extract clean DNA reads or k-mers
Use the provided Python scripts to convert the FASTQ data into formats compatible with your PostgreSQL extension.

*Generate DNA sequences:*
- ```python extract_dna.py data/SRR026760.fastq 0 > dna_36.txt```

*Generate Q-Kmers:*
- ```python generate_qkmer.py data/SRR026760.fastq 32 > qkmeroutput.txt```

## Summary
Using synthetic sequences ensures your code behaves correctly in controlled cases, while real-world SRA datasets confirm that your extension scales to real genomic data and handles biological noise. Together, these tests provide complete coverage for validating the k-mer/q-mer operators, SP-GiST indexing, and overall functionality of the DNA extension.
