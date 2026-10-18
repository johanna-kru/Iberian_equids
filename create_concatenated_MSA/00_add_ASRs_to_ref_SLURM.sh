#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --job-name=ASR2ref
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu 18000
#SBATCH -c 3
#SBATCH --partition compute


###	ALIGN ASRs (ANCIENT SEQUENCES RECONSTRUCTIONS) TO MULTIPLE SEQUENCE ALIGNMENT OF REFERENCE SEQUENCES


######
### Arguments
######

refalns="$1"
asrs="$2"
outdir="$3"
gois="$4"
samples="$5"
goi=$(awk -v num=$SLURM_ARRAY_TASK_ID 'NR==num' ${gois})

# there should be an existing MSA for every goi, otherwise just create a file with a UniProt reference


#######
### Modules
######

# module load seqtk/1.3-GCC-11.2.0

########
### Read in files
######

echo "#################### Reading in files for ${goi} ##############################"


		echo "### Parsing files for ${goi} ###"


######
### Create files to align
######

		mkdir -p "$outdir/tmp_${goi}"


		### concatenate sequences to align

		touch "$outdir/tmp_${goi}/${goi}_unaligned.prefasta" 	#{goi}_unaligned in {} otherwise file not correctly written, can be used but not visible


		cat "$asrs/"*${goi}*".fasta" >> "$outdir/tmp_${goi}/${goi}_unaligned.prefasta"
		echo "$outdir/tmp_${goi}/${goi}_unaligned.prefasta created."

		# exchange gaps for X for better alignment
		sed '/^>/! s/-/X/g' "$outdir/tmp_${goi}/${goi}_unaligned.prefasta" > "$outdir/tmp_${goi}/${goi}_unaligned.fasta"
		
	



		### selected individuals for reference dataset
		cat "${refalns}/"*${goi}*".fasta" | grep ^">" | grep -f "$samples" | cut -d">" -f2 > "$outdir/tmp_${goi}/${goi}_headers.txt"

		headers="$outdir/tmp_${goi}/${goi}_headers.txt"
 
  		awk -v var="$headers" 'BEGIN{while((getline<var)>0)l[">"$1]=1}/^>/{f=l[$1]}f' "${refalns}/"*${goi}*".fasta" > "$outdir/tmp_${goi}/${goi}_selected_ref.fasta"
		
		msafile=$(echo "$outdir/tmp_${goi}/${goi}_selected_ref.fasta")

#		seqtk subseq "$outdir/tmp_${goi}/${goi}_unaligned.prefasta" "$outdir/tmp_${goi}/${goi}_headers.txt" > "$outdir/tmp_${goi}/${goi}_unaligned.fasta"
# seqtk cannot deal with long headers!



### align inmsa because mafft complains they are not aligned

	echo "align inmsa for ${goi} because mafft complains they are not aligned" >&2
		mafft \
			--maxiterate 1000 \
			"$msafile" \
			> "$outdir/tmp_${goi}/${goi}_inmsa.fasta"
# older versions of mafft that does not have addfragemnts and complains about shell config 

#######
### Align sequences to ref
#######

	echo "Align sequences to ref for ${goi}" >&2

		maffadd="$outdir/tmp_${goi}/${goi}_unaligned.fasta"
		maffin="$outdir/tmp_${goi}/${goi}_inmsa.fasta"

		mafft --thread 3 --maxiterate 1000 --add $maffadd $maffin > "$outdir/${goi}_addfrag.fasta"

			rm -r "$outdir/tmp_${goi}"
