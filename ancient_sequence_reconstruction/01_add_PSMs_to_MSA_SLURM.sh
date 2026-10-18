#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --job-name=add_PSMs
#SBATCH --mem-per-cpu 10000
#SBATCH -c 3
#SBATCH --partition compute


######
### Arguments
######

refalns=$1
outdir=$2
gois=$3
samples=$4
sample=$(awk -v num=$SLURM_ARRAY_TASK_ID 'NR==num' $samples)
addalns=$5

# there should be an existing MSA for every goi, otherwise just create a file with a UniProt reference

#######
### Modules
######

# module load seqtk/1.3-GCC-11.2.0


########
### Read in files
######

echo "Reading in files for $sample"

	while read goi; do 

		echo "Parsing files for $goi"

		msafile=$(echo "$refalns/"*${goi}*".fasta")

		addfile1=$(echo "$addalns/"*${sample}"_"*".fasta")

			if [ "$6" != "" ]; then
				addfile2=$(echo "$6/"*${goi}*".fasta") 
			else
				addfile2=""
			fi


			if [ "$7" != "" ]; then
				addfile3=$(echo "$7/"*${goi}*".fasta") 
			else
				addfile3=""
			fi


			if [ "$8" != "" ]; then
				addfile4=$(echo "$8/"*${goi}*".fasta") 
			else
				addfile4=""
			fi




######
### Create files to align
######

		mkdir -p "$outdir/$sample/tmp_$goi"


		### concatenate sequences to align

		touch "$outdir/$sample/tmp_$goi/${goi}_unaligned.prefasta"

		cat $addfile1 $addfile2 $addfile3 $addfile4 >> "$outdir/$sample/tmp_$goi/${goi}_unaligned.prefasta"
		echo "$outdir/$sample/tmp_$goi/${goi}_unaligned.prefasta created."

		#{goi}_unaligned in {} otherwise file not correctly written, can be used but not visible

		### extract only goi
		cat "$outdir/$sample/tmp_$goi/${goi}_unaligned.prefasta" | grep ^">" | grep $goi | cut -d">" -f2 > "$outdir/$sample/tmp_$goi/${goi}_headers.txt"

		headers="$outdir/$sample/tmp_$goi/${goi}_headers.txt"
 
  		awk -v var="$headers" 'BEGIN{while((getline<var)>0)l[">"$1]=1}/^>/{f=l[$1]}f' "$outdir/$sample/tmp_$goi/${goi}_unaligned.prefasta" > "$outdir/$sample/tmp_$goi/${goi}_unaligned.fasta"

# seqtk subseq "$outdir/$sample/tmp_$goi/${goi}_unaligned.prefasta" "$outdir/$sample/tmp_$goi/${goi}_headers.txt" > "$outdir/$sample/tmp_$goi/${goi}_unaligned.fasta"
# seqtk cannot deal with long headers



### align inmsa because mafft complains they are not aligned
	echo MAFFT_BINARIES $MAFFT_BINARIES
	unset MAFFT_BINARIES
	ls -l $msafile

	echo "align inmsa for $goi because mafft complains they are not aligned" >&2
		~/miniconda3/bin/mafft \
			--maxiterate 1000 \
			"$msafile" \
			> "$outdir/$sample/tmp_$goi/${goi}_inmsa.fasta"

#######
### Align sequences to ref
#######

	echo "Align sequences to ref for $goi" >&2

		maffadd="$outdir/$sample/tmp_$goi/${goi}_unaligned.fasta"
		maffin="$outdir/$sample/tmp_$goi/${goi}_inmsa.fasta"
		
		echo MAFFT_BINARIES $MAFFT_BINARIES
		unset MAFFT_BINARIES
		~/miniconda3/bin/mafft --thread 3 --maxiterate 1000 --addfragments $maffadd $maffin > "$outdir/$sample/${goi}_addfrag.fasta"

		# 	--keeplength \
		# addfragments does not exist in this version (v7.490)!
		# conda version (v7.526) has addfragments
		# --op gap open penalty, default 1.53
		# for --add use --op 3 and --localpair
		rm -r "$outdir/$sample/tmp_$goi"

	done < "$gois"