#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --job-name=IQtree
#SBATCH --mem-per-cpu 14000
#SBATCH -c 1
#SBATCH --partition compute

######
### Arguments
######

indir=$1
outdir=$2
samples=$3

sample=$(awk -v num=$SLURM_ARRAY_TASK_ID 'NR==num' $samples)
inaln="$indir/$sample"
samplestem=$(echo "$sample" | cut -d"." -f1)
outstem="IQtree_$samplestem"

#######
### Load modules
#######


module load bio/iqtree/2.3.6




######
### Prepare modules and directory
######


    sampledir="$outdir/$samplestem"

    mkdir -p "$sampledir"

    cd "$sampledir"



######
### Calculate Maximum Likelihood tree with IQtree
######

trd=1

iqtree2 \
    -s "${inaln}" \
    -m MFP \
    -mset LG+F+G,WAG+F+G,JTT+F+G,GTR20 \
    -madd LG+C20+F+G,LG+C10+F+G,LG+C30+F+G,LG+C40+F+G,LG+C50+F+G,LG+C60+F+G,C10,C20,C30,C40,C50,C60,EX2,EX3,EHO,LG4M,LG4X \
    -nt $trd \
    -nstop 500 \
    -bb 1000 \
    -alrt 5000 \
    -pre "${sampledir}/${outstem}"









    # /scratch/devel/jkrueger/bin/iqtree-1.6.12-Linux/bin/iqtree \
    #         -s $inaln  \
    #         -spp $partfile \
    #         -alrt 5000 \
    #         -bb 5000 \
	# 	    -nstop 500 \
	#         -nm 10000 \
    #         -wsplits \
    #         -nt $trd \
    #         -pre $outdir/tmp_$SLURM_ARRAY_TASK_ID/$outstem


