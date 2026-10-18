#!/usr/bin/env python3

### created with help of a large language 

import sys
from Bio import SeqIO

def main(fasta_file, header_file, output_file):
    # Read headers from header list
    with open(header_file) as f:
        ordered_headers = [line.strip() for line in f if line.strip()]

    # Parse the FASTA and build a dictionary {header: record}
    records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Write matched records in the order of the header list
    with open(output_file, "w") as out:
        for header in ordered_headers:
            if header in records:
                SeqIO.write(records[header], out, "fasta")
            else:
                sys.stderr.write(f"Warning: Header '{header}' not found in FASTA file\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python sort_fasta_by_headers.py input.fasta header_list.txt output.fasta")
        sys.exit(1)

    fasta_file     = sys.argv[1]
    header_file    = sys.argv[2]
    output_file    = sys.argv[3]

    main(fasta_file, header_file, output_file)
