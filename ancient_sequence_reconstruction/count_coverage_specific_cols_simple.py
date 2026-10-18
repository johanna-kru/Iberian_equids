### written with support of large language models

from collections import Counter

def read_fasta(filepath):
    sequences = {}
    with open(filepath, "r") as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences[header] = "".join(seq_lines)
                header = line[1:]  # remove ">"
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            sequences[header] = "".join(seq_lines)
    return sequences


def column_composition(msa_file, column_number):
    sequences = read_fasta(msa_file)

    # Filter sequences containing "ricScore"
    filtered = {
        h: s for h, s in sequences.items()
        if "ricScore" in h
    }

    if not filtered:
        print("No sequences found with 'ricScore' in header.")
        return

    # Convert to 0-based index
    col_idx = column_number - 1

    residues = []
    for header, seq in filtered.items():
        if col_idx >= len(seq):
            print(f"Warning: Sequence '{header}' is shorter than column {column_number}")
            continue
        residues.append(seq[col_idx])

    counts = Counter(residues)

    print(f"\nColumn {column_number} composition (absolute counts):")
    for residue, count in sorted(counts.items()):
        print(f"{residue}: {count}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Get column composition from MSA")
    parser.add_argument("msa", help="MSA file in FASTA format")
    parser.add_argument("column", type=int, help="Column number (1-based)")

    args = parser.parse_args()

    column_composition(args.msa, args.column)
