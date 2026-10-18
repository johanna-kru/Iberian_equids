#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --job-name=getcoverage
#SBATCH --mem-per-cpu 8000
#SBATCH -c 1
#SBATCH --partition=compute

######
### Arguments
######

indir=$1
outdir=$2
gois=$3
samples=$4
sample=$(awk -v num=$SLURM_ARRAY_TASK_ID 'NR==num' $samples)


######
### Prepare calling consensus
######

	mkdir -p $outdir


########
### Read in files
######


	while read goi; do 

		echo "Calculate coverage of protein $goi for $sample"

			inmsa=$indir/$sample/*${goi}*.fasta

			python calculate_coverage.py $inmsa > $outdir/"coverage_${goi}_${sample}.tsv"


	done < "$gois"