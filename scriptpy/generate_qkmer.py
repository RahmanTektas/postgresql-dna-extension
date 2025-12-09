import sys
import random


COMPATIBLE_IUPAC = {
    'A': ['R', 'M', 'W', 'V', 'H', 'D'],  # Ambiguities involving A
    'C': ['Y', 'M', 'S', 'V', 'H', 'B'],  # Ambiguities involving C
    'G': ['R', 'S', 'K', 'V', 'D', 'B'],  # Ambiguities involving G
    'T': ['Y', 'W', 'K', 'H', 'D', 'B'],  # Ambiguities involving T
}


def dna_to_qkmer(kmer, prob_ambig=0.3, prob_N=0.05):
    """
    Convert a DNA k-mer (A/C/G/T) into a q-mer by introducing biologically valid
        IUPAC ambiguity symbols.

    Ambiguity substitution respects compatibility rules:
        - Each ambiguous symbol is only used if the original base belongs to
        its definition. For example:
            Y (C/T) may only replace C or T
            R (A/G) may only replace A or G
            K (G/T) may only replace G or T
            etc.
        - Therefore, all generated ambiguity codes reflect valid nucleotide sets.

    Replacement logic:
        • With probability `prob_N`, the base is replaced with 'N'
        (fully ambiguous / unknown).
        • Otherwise, with probability `prob_ambig`, the base is replaced with
        a random compatible ambiguity symbol among:
            R, Y, M, W, S, K, V, H, D, B
        • If no ambiguity is introduced, the original base is kept unchanged.

    Args:
        kmer (str): Input DNA k-mer consisting of A, C, G, T.
        prob_ambig (float): Probability of replacing a base with a compatible
                            IUPAC ambiguity code.
        prob_N (float): Probability of forcing the base to N.

    Returns:
        str: The generated q-mer containing valid DNA/IUPAC characters.
    """
    q = []

    for b in kmer.upper():
        r = random.random()

        # Force N
        if r < prob_N:
            q.append('N')
            continue

        # Introduce a valid ambiguity
        if r < prob_N + prob_ambig and b in COMPATIBLE_IUPAC:
            q.append(random.choice(COMPATIBLE_IUPAC[b]))
        else:
            q.append(b)

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
        print("Usage: python generate_qkmer.py <fastq> <k> [max_reads]", file=sys.stderr)
        sys.exit(1)

    fastq = sys.argv[1]
    k = int(sys.argv[2])
    max_reads = int(sys.argv[3]) if len(sys.argv) >= 4 else None

    fastq_to_qkmer(fastq, k, max_reads)
