import sys

def generate_dna(seq, k):
    # on garde la signature, mais on ignore k
    seq = seq.strip().upper()
    # si tu veux ignorer les reads avec N, on garde ce comportement
    if seq and "N" not in seq:
        yield seq

def fastq_to_dna(fastq_path, k, max_reads=None):
    # k ne sert plus vraiment mais on le garde pour la compatibilité
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

            # une ligne = un ADN (un read complet)
            for dna in generate_dna(seq, k):
                print(dna)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_dna.py <fastq> <k> [max_reads]", file=sys.stderr)
        sys.exit(1)

    fastq = sys.argv[1]
    k = int(sys.argv[2])
    max_reads = int(sys.argv[3]) if len(sys.argv) >= 4 else None
    fastq_to_dna(fastq, k, max_reads)
