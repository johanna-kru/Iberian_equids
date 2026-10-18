### written with support of large language models
from Bio import AlignIO
from collections import Counter
import sys

def consensus(mult_seq_alignment):
    consensus_seq = ""
    length = mult_seq_alignment.get_alignment_length()

    for i in range(length):
        # Extract the bases at the current position in the alignment
        bases = mult_seq_alignment[:, i]
        # Count the frequency of each base
        base_counts = Counter(bases)
        # Ignore gaps in the count
        del base_counts["-"]
        # If all were gaps, add a gap to the consensus sequence
        if len(base_counts) == 0:
            consensus_seq += "-"
        else:
            # Find the base with the maximum count
            max_base = max(base_counts, key=base_counts.get)
            consensus_seq += max_base

    return consensus_seq

# Read the multiple sequence alignment file
alignment_file = sys.argv[1]
alignment = AlignIO.read(alignment_file, "fasta")

# Compute the consensus sequence
consensus_sequence = consensus(alignment)

# Print the consensus sequence
print(">consensus\n" + consensus_sequence)
