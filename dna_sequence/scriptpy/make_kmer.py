import sys

def generate_kmers(seq, k):
    seq = seq.strip().upper()
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i:i+k]
        # on ignore les séquences avec N (base inconnue)
        if "N" not in kmer:
            yield kmer

def fastq_to_kmers(fastq_path, k, max_reads=None):
    k = int(k)
    count_reads = 0
    with open(fastq_path, "r") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if not qual:
                break

            count_reads += 1
            if max_reads is not None and count_reads > max_reads:
                break

            for kmer in generate_kmers(seq, k):
                # une ligne = un kmer
                print(kmer)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python make_kmers.py <fastq> <k> [max_reads]", file=sys.stderr)
        sys.exit(1)

    fastq = sys.argv[1]
    k = int(sys.argv[2])
    max_reads = int(sys.argv[3]) if len(sys.argv) >= 4 else None
    fastq_to_kmers(fastq, k, max_reads)
