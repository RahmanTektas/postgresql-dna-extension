import sys
import random

def dna_to_qkmer(kmer, prob_ambig=0.3, prob_N=0.05):
    """
    Transforme un k-mer (A/C/G/T) en qkmer (A,C,G,T,R,Y,N seulement).

    - prob_ambig = proba qu'une base devienne R ou Y
    - prob_N     = proba qu'une base devienne N (joker total)
    """
    q = []
    for b in kmer:
        b = b.upper()
        r = random.random()

        # De temps en temps, on met un N (wildcard total)
        if r < prob_N:
            q.append('N')
            continue

        # Sinon on peut mettre une ambiguïté R ou Y
        if b in ('A', 'G'):
            # purine -> soit A/G soit R
            if r < prob_N + prob_ambig:
                q.append('R')
            else:
                q.append(b)
        elif b in ('C', 'T'):
            # pyrimidine -> soit C/T soit Y
            if r < prob_N + prob_ambig:
                q.append('Y')
            else:
                q.append(b)
        else:
            # au cas où (N ou autre), on met N
            q.append('N')

    return "".join(q)

def generate_dna(seq, k):
    seq = seq.strip().upper()
    for i in range(0, len(seq) - k + 1):
        dna = seq[i:i+k]
        # on ignore les séquences avec N (base inconnue dans le read)
        if "N" not in dna:
            yield dna

def fastq_to_qkmer(fastq_path, k, max_reads=None):
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

            for dna in generate_dna(seq, k):
                qkmer = dna_to_qkmer(dna)
                print(qkmer)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_qkmer.py <fastq> <k> [max_reads]", file=sys.stderr)
        sys.exit(1)

    fastq = sys.argv[1]
    k = int(sys.argv[2])
    max_reads = int(sys.argv[3]) if len(sys.argv) >= 4 else None
    fastq_to_qkmer(fastq, k, max_reads)
