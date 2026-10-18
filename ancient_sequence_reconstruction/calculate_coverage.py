### written with support of large language models

from Bio import SeqIO
import sys

def calculate_coverage(fasta_file):
    alignment = SeqIO.parse(fasta_file, "fasta")
    coverage = [0]*len(next(alignment))
    for record in alignment:
        for i, aa in enumerate(record.seq):
            if aa != '-':
                coverage[i] += 1
    return coverage

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    coverage = calculate_coverage(fasta_file)
    print("position\tcoverage")
    for i, cov in enumerate(coverage, start=1):
        print(f"{i}\t{cov}")
