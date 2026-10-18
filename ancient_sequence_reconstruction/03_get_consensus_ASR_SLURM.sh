#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --job-name=consensus_ASR
#SBATCH --mem-per-cpu 8000
#SBATCH -c 1
#SBATCH --partition compute

### Make sure reference sequcnes are deleted from input alignments. Use selection script, e.g. the ones in general

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

		echo "Create consensus sequence of protein $goi for $sample"

			inmsa=$indir/$sample/*${goi}*.fasta

			python get_consensus.py $inmsa > $outdir/"consensus_${goi}_${sample}.prefasta"

			cat $outdir/"consensus_${goi}_${sample}.prefasta" | sed "s/>/>${goi}_${sample}_/g" > $outdir/"consensus_${goi}_${sample}.fasta"

			rm $outdir/"consensus_${goi}_${sample}.prefasta"

	done < "$gois"