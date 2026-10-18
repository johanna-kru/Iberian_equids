#!/bin/bash

#SBATCH --time=00:05:00
#SBATCH --job-name=count_aa
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu 8000
#SBATCH -c 1
#SBATCH --qos=test
#SBATCH --partition=genD



######
### Arguments
######

indir="consensus"
outdir="consensus"


######
### Prepapre directories
######


mkdir -p $outdir

ls $indir/*.fasta | rev | cut -d"/" -f1 | rev > $indir/"filenames.txt"

filenames=$indir/"filenames.txt"


######
### Count amino acids
######

rm $outdir/"aa_counts.tsv"
touch $outdir/"aa_counts.tsv"

    while read infile; do
        
        sample=$(echo $indir/$infile | rev | cut -d"/" -f1 | rev | cut -d"_" -f3 | cut -d"." -f1 )
        echo $sample

        goi=$(echo $indir/$infile | rev | cut -d"/" -f1 | rev |  cut -d"_" -f2)
        echo $goi

        aanum=$(cat $indir/$infile | grep -v ^">" | grep -o [A-Z] | wc -l)

        echo -e "$sample\t$goi\t$aanum" >> $outdir/"aa_counts.tsv"


    done < "$filenames"