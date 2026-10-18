### written with support of large language models

## for coverage of interesting single amino acis variants (SAVs)

import os
import re
import csv
from collections import Counter, defaultdict


# ---------- FASTA READER ----------
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
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            sequences[header] = "".join(seq_lines)
    return sequences


# ---------- GROUP SAMPLES ----------
def group_samples(sample_list):
    """
    Groups samples:
    - P1327A + P1327B -> P1327AB
    - Special: P2855A, P2856A, P2857A -> P2855A-2857A
    """

    grouped = defaultdict(list)

    for sample in sample_list:
        match = re.match(r"(P\d+)([A-Z])", sample)
        if not match:
            continue
        prefix, letter = match.groups()
        grouped[prefix].append(letter)

    final_groups = {}

    for prefix, letters in grouped.items():
        letters = sorted(set(letters))  # avoid duplicates

        # ---- SPECIAL CASE ----
        if prefix in ["P2855", "P2856", "P2857"]:
            key = "P2855A-2857A"
            final_groups.setdefault(key, [])
            for l in letters:
                final_groups[key].append(f"{prefix}{l}")
            continue

        # ---- NORMAL CASE ----
        if len(letters) == 1:
            key = f"{prefix}{letters[0]}"
        else:
            key = f"{prefix}{''.join(letters)}"

        final_groups[key] = [f"{prefix}{l}" for l in letters]

    return final_groups


# ---------- COLUMN COUNT ----------
def get_column_counts(fasta_file, column):
    if not os.path.exists(fasta_file):
        return None

    seqs = read_fasta(fasta_file)

    # Filter sequences containing "ricScore"
    seqs = {h: s for h, s in seqs.items() if "ricScore" in h}

    if not seqs:
        return Counter()

    col_idx = column - 1
    residues = []

    for seq in seqs.values():
        if col_idx < len(seq):
            aa = seq[col_idx]
            if aa != "-":  # ignore gaps
                residues.append(aa)

    return Counter(residues)


# ---------- FORMAT COUNTS ----------
def format_counts(counter):
    if not counter or sum(counter.values()) == 0:
        return "NA"
    return ";".join([f"{aa}:{count}" for aa, count in sorted(counter.items())])


# ---------- MAIN ----------
def main(sample_file, protein_csv, output_file):
    base_path = "aln_to_ref_addfragments"

    # Read sample list
    with open(sample_file) as f:
        samples = [line.strip() for line in f if line.strip()]

    # Group samples automatically
    grouped_samples = group_samples(samples)

    # Read protein-position CSV
    with open(protein_csv) as f:
        reader = csv.DictReader(f)
        protein_positions = list(reader)

    # Prepare output header
    header = ["protein", "position"] + list(grouped_samples.keys())

    with open(output_file, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(header)

        for entry in protein_positions:
            protein = entry["protein"]
            position = int(entry["position"])

            row = [protein, position]

            for group_name, sub_samples in grouped_samples.items():
                total_counter = Counter()

                for sub in sub_samples:
                    fasta_path = os.path.join(
                        base_path,
                        sub,
                        f"{protein}_addfrag.fasta"
                    )

                    counts = get_column_counts(fasta_path, position)

                    if counts:
                        total_counter.update(counts)

                row.append(format_counts(total_counter))

            writer.writerow(row)


# ---------- RUN ----------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate residue composition table from MSAs"
    )
    parser.add_argument("--samples", required=True, help="Sample list (txt)")
    parser.add_argument("--proteins", required=True, help="Protein-position CSV")
    parser.add_argument("--output", required=True, help="Output CSV file")

    args = parser.parse_args()

    main(args.samples, args.proteins, args.output)