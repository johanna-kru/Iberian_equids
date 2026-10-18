#!/bin/bash

### EXAMPLE PARAMETER FILE. THESE FILES WERE OFTEN USED TO SUBMIT JOBS TO SCRIPTS NAMED "somename_SLURM.sh". 

######
### PARAMETER FILE 
######
refalns="ref_MSAs/"
outdir="aln_to_ref_addfragments/"
gois="gois.txt"
samples="samples.txt"
addalns="accepted_peptides/"



no_arrays="1-23"
# if only one "1" or "5" (line number of item in list), if many "1-12"



	sbatch \
		--array=$no_arrays \
		--output=log/array_%A-%a.out \
		--error=log/array_%A-%a.err \
		add_fragments_to_MSA_SLURM.sh \
		$refalns \
		$outdir \
		$gois \
		$samples \
		$addalns
