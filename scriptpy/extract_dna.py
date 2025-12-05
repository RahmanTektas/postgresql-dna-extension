import sys

def generate_dna(seq, k):
    """
    Yield a cleaned DNA sequence if valid.

    This function is intentionally simple: it returns the full read
    only if:
        - the sequence is non-empty
        - it contains no 'N' (ambiguous nucleotide)

    Note:
        The parameter `k` is accepted for interface consistency with
        other extractors, but not used here since we return the whole
        read as-is.

    Args:
        seq (str): The sequence read from FASTQ.
        k (int): k-mer length (unused).
    
    Yields:
        str: The cleaned DNA sequence.
    """
    seq = seq.strip().upper()

    # Only yield if the read is fully defined (no ambiguous bases)
    if seq and "N" not in seq:
        yield seq


def fastq_to_dna(fastq_path, k, max_reads=None):
    """
    Extract full-length DNA sequences from a FASTQ file.

    FASTQ files are read 4 lines at a time:
        1. @header
        2. sequence
        3. +
        4. quality scores

    This function prints each cleaned DNA read directly to stdout.

    Args:
        fastq_path (str): Path to the FASTQ input file.
        k (int): Placeholder k-mer length (not used in this version).
        max_reads (int or None): Optional limit on number of reads processed.
    """
    k = int(k)  # kept for consistency with qkmer/kmer tools
    count_reads = 0

    with open(fastq_path, "r") as f:
        while True:
            # Read one FASTQ record
            header = f.readline()
            if not header:
                break  # End of file reached

            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            if not qual:
                break  # Incomplete FASTQ record

            count_reads += 1

            # Stop early if a read limit was set
            if max_reads is not None and count_reads > max_reads:
                break

            # Emit the full sequence if valid
            for dna in generate_dna(seq, k):
                print(dna)


if __name__ == "__main__":
    # Expected usage:
    #   python extract_dna.py <fastq> <k> [max_reads]
    if len(sys.argv) < 3:
        print("Usage: python extract_dna.py <fastq> <k> [max_reads]", file=sys.stderr)
        sys.exit(1)

    fastq = sys.argv[1]
    k = int(sys.argv[2])
    max_reads = int(sys.argv[3]) if len(sys.argv) >= 4 else None

    fastq_to_dna(fastq, k, max_reads)
