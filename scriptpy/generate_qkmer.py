import sys
import random

def dna_to_qkmer(kmer, prob_ambig=0.3, prob_N=0.05):
    """
    Convert a DNA k-mer into a q-mer (ambiguous nucleotide representation).

    For each nucleotide:
      - With probability prob_N: replace with 'N' (fully ambiguous).
      - Otherwise:
          - If base is A or G (purines):
                With prob prob_ambig: replace with 'R'
                Else: keep original base
          - If base is C or T (pyrimidines):
                With prob prob_ambig: replace with 'Y'
                Else: keep original base
          - Any unexpected base defaults to 'N'

    Args:
        kmer (str): Input DNA k-mer (A, C, T, G).
        prob_ambig (float): Probability of introducing IUPAC ambiguity (R/Y).
        prob_N (float): Probability of forcing the base to 'N'.

    Returns:
        str: The generated q-kmer.
    """
    q = []
    for b in kmer:
        b = b.upper()
        r = random.random()

        # Force N with small probability
        if r < prob_N:
            q.append('N')
            continue

        # Purines (A/G)
        if b in ('A', 'G'):
            if r < prob_N + prob_ambig:
                q.append('R')  # Ambiguous purine
            else:
                q.append(b)

        # Pyrimidines (C/T)
        elif b in ('C', 'T'):
            if r < prob_N + prob_ambig:
                q.append('Y')  # Ambiguous pyrimidine
            else:
                q.append(b)

        # Any unknown character → N
        else:
            q.append('N')

    return "".join(q)

def generate_dna(seq, k):
    """
    Generate all clean k-mers (no 'N') from a DNA sequence.

    Args:
        seq (str): Raw DNA sequence (may contain lowercase).
        k (int): Length of k-mers to extract.

    Yields:
        str: Clean DNA k-mer (only A, C, G, T).
    """
    seq = seq.strip().upper()

    # Extract sliding window k-mers
    for i in range(0, len(seq) - k + 1):
        dna = seq[i:i+k]

        # Skip k-mers containing ambiguous bases already
        if "N" not in dna:
            yield dna

def fastq_to_qkmer(fastq_path, k, max_reads=None):
    """
    Convert all valid DNA substrings from a FASTQ file into q-kmers.

    FASTQ reading is done 4 lines at a time:
        1. @header
        2. sequence
        3. +
        4. quality scores

    Args:
        fastq_path (str): Path to the FASTQ file.
        k (int): k-mer length.
        max_reads (int or None): If set, stop after processing this many reads.

    Output:
        Prints each q-kmer to stdout.
    """
    k = int(k)
    count_reads = 0

    with open(fastq_path, "r") as f:
        while True:
            # Read one FASTQ record
            header = f.readline()
            if not header:
                break  # End of file

            seq = f.readline()
            plus = f.readline()
            qual = f.readline()

            if not qual:
                break  # Malformed FASTQ or unexpected end

            count_reads += 1

            # Optional limiting of processed reads
            if max_reads is not None and count_reads > max_reads:
                break

            # Extract clean k-mers and convert them
            for dna in generate_dna(seq, k):
                qkmer = dna_to_qkmer(dna)
                print(qkmer)

if __name__ == "__main__":
    # Expect: python extract_qkmer.py <fastq> <k> [max_reads]
    if len(sys.argv) < 3:
        print("Usage: python extract_qkmer.py <fastq> <k> [max_reads]", file=sys.stderr)
        sys.exit(1)

    fastq = sys.argv[1]
    k = int(sys.argv[2])
    max_reads = int(sys.argv[3]) if len(sys.argv) >= 4 else None

    fastq_to_qkmer(fastq, k, max_reads)
