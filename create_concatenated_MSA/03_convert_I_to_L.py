#use with PYTHON/3.8.6

from Bio import Seq, SeqIO
import sys

if len(sys.argv) != 3:
    print("Usage:", sys.argv[0], "<input_fasta_file> <asr_list_file>")
    exit(1)

input_filename = sys.argv[1]
asr_filename = sys.argv[2]

# Read the ASR file
asr_names = []
with open(asr_filename) as asr_file:
    for line in asr_file:
        line = line.strip()
        if line.startswith("#"):
            continue # commented out line
        if line.startswith(">"):
           line = line[1:] # remove preceding '>' symbols
        if line:
            asr_names.append(line)

# Read the fasta file and split the records between ASR and reference
asr_records = []
ref_records = []
fasta_records = SeqIO.parse(input_filename, 'fasta')
for record in fasta_records:
    name = record.id
    sequence = str(record.seq)
    #print(name)
    if name in asr_names:
        asr_records.append(record)
    else:
        ref_records.append(record)

# Error handling: some of the ASR names in the input ASR file were not found in the fasta file
if len(asr_records) != len(asr_names):
    print("Some ASRs were not found in the fasta file:")
    for name in asr_names:
        found = False
        for record in asr_records:
            if name == record.id:
                found = True
                break
        if not found:
            print(" *", name)
    exit(2)

# Read non-empty positions from ASR records
relevant_positions = set()
for record in asr_records:
    for pos, char in enumerate(record.seq):
        if char == 'I' or char == 'L' or char == 'i' or char == 'l':
           relevant_positions.add(pos)
print("Found", len(relevant_positions), "relevant positions in the ASR sequences")



# Change all 'I' to 'L' in relevant positions from ASR records
from Bio.Seq import Seq

for record in asr_records:
    sequence = list(str(record.seq))
    for p in relevant_positions:
        if p > len(sequence):
            break
        if sequence[p] == 'I' or sequence[p] == 'i':
              sequence[p] = 'L'
    record.seq = Seq(''.join(sequence))



# Change all 'I' to 'L' in relevant positions from references records
for record in ref_records:
    sequence = list(str(record.seq))
    for p in relevant_positions:
        if p > len(sequence):
            break
        if sequence[p] == 'I' or sequence[p] == 'i':
              sequence[p] = 'L'
    record.seq = Seq(''.join(sequence))



# Write out the modified fasta file
input_filename_minus_extension = input_filename.rsplit('.', 1)[0]
output_filename = "{}_replacedI.fasta".format(input_filename_minus_extension)
SeqIO.write(asr_records + ref_records, output_filename, "fasta")

