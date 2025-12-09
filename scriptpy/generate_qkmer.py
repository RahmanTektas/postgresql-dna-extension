import sys
import random


COMPATIBLE_IUPAC = {
    'A': ['R', 'M', 'W', 'V', 'H', 'D'],  # Ambiguities involving A
    'C': ['Y', 'M', 'S', 'V', 'H', 'B'],  # Ambiguities involving C
    'G': ['R', 'S', 'K', 'V', 'D', 'B'],  # Ambiguities involving G
    'T': ['Y', 'W', 'K', 'H', 'D', 'B'],  # Ambiguities involving T
}


def dna_to_qkmer(kmer, prob_ambig=0.3, prob_N=0.05):
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
    seq = seq.strip().upper()

    # Extract sliding window k-mers
    for i in range(0, len(seq) - k + 1):
        dna = seq[i:i+k]

        # Skip k-mers containing ambiguous bases already
        if "N" not in dna:
            yield dna

def fastq_to_qkmer(fastq_path, k, max_reads=None):
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
