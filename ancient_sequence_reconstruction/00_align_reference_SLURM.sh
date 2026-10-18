#!/bin/bash

#SBATCH --time=1-23:59:00
#SBATCH --job-name=aln_predictions
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu 30000
#SBATCH -c 3
#SBATCH --qos=vlong
#SBATCH --partition=genD

#####################################################################

######
### Modules
######



######
### Arguments
######

indir=$1 
outdir=$2
gois=$3
prefix=$4

ntrd=3

######
### Prepare run
######


mkdir -p $outdir

goi=$(awk -v num=$SLURM_ARRAY_TASK_ID 'NR==num' $gois)

######
### Create alignment
######

cat \
$indir/*/"proteins"/*/$goi/*"_private_CDS_Prot_sequences.fasta" \
"/ref_proteins/${goi}_reference.fasta" \
"/CerSimSim1.0/protein_predictions/${goi}/Cersimsim1_private_CDS_Prot_sequences_tsv.fasta" \
| mafft --thread $ntrd --unalignlevel 0.0 --maxiterate 100 --amino --reorder --globalpair --bl 80 - > $outdir/$prefix"_"$goi"_Prot.aln.fasta"

# globalpair: All pairwise alignments are computed with the Needleman-Wunsch algorithm. More accurate but slower than --6merpair. Suitable for a set of globally alignable sequences. Applicable to up to ~200 sequences. A combination with --maxiterate 1000 is recommended (G-INS-i). Default: off (6mer distance is used)
  # maxiterate: number cycles of iterative refinement are performed. Default: 0
  # bl: BLOSUM Matix, 80 because of conserved sequences
  # --reorder: "Output order: aligned. Default: off (inputorder)""
  # unalignlevel 0 tries to align X, the higher the easier gaps formed around non-homologous sequences

