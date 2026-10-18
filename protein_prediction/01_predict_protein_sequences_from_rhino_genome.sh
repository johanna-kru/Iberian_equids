#!/bin/bash

# Code used to translate sequences directly from the reference genoms. Any reference genome can be queried using this code by using the pertinent gtf file for the given assembly.


######
### Modules
######

module load GCC/7.3.0-2.30
module load XZ/5.2.4-GCCcore-7.3.0
module load BCFtools/1.14-GCC-11.2.0
module load SAMtools/1.14-GCC-11.2.0
module load zlib/1.2.8-GCCcore-5.4.0
module load Python/3.8.6-GCCcore-10.2.0



######
### Arguments
######



work_dir=""
pydir="directory_with_Combine_CDS_concatenator.py"
refgen=${work_dir}"/GCF_000283155.1_CerSimSim1.0_genomic.fna"
gois="cersimsim1.txt"





######
### Translating genomic data
######

### loop over gois

	while read goi; do
	outdir="outdir/${goi}"
	mkdir -p $outdir

		GFF=${work_dir}"/${goi}_CerSimSim1.0.refseq.gtf"

		# This line now is just retrieving one gene. You may use a list of gene identifiers of your given organism. It is easier if you just build a list in a separate file, and then feed it to this program one by one by using a loop.
		Tupaia_names=$(cat $GFF | awk '{print $12}' | sed -e 's/\"//g' -e 's/\;//g' | sort | uniq)
		#Tupaia_names="XM_014585727.1"


		echo "Building tsv file:"

		for identifier in ${Tupaia_names[@]}; do

			echo "for "$identifier
			# chunk cortado, solo salvando: grep -v start_codon | grep -v stop_codon |sort -nk4 
			cat $GFF | grep $identifier | grep CDS | grep -v start_codon | grep -v stop_codon | sort -nk4 | while read line; do 

				scaffold=$(echo $line | awk '{print $1}');
				start_co=$(echo $line| awk '{print $4}');
				stop_co=$(echo $line | awk '{print $5}');
				orientation=$(echo $line | awk '{print $7}');
				gene=$(echo $line | awk '{print $10}' | sed -e 's/\;//g' -e 's/\"//g');
				tag=$(echo $line | awk '{print $12}' | sed -e 's/\;//g' -e 's/\"//g');

				samtools faidx $refgen ${scaffold}:${start_co}-${stop_co} | \
				/scratch_isilon/groups/compgen/jkrueger/bin/bioawk-master/bioawk  -cfastx \
					-vtag=$tag \
					-vstart_co=$start_co \
					-vstop_co=$stop_co \
					-vscaffold=$scaffold \
					-vorientation=$orientation \
					-vgene=$gene \
				'{OFS="\t";print scaffold, start_co-1, stop_co, orientation, tag, $seq, gene}' | \
				sort -k1,1 -k2,2n >> ${outdir}/"Cersimsim1.private_exon_sequences.tsv";
				
			done;
		done;

		python3 ${pydir}/Combine_CDS_concatenator.py ${outdir}/"Cersimsim1.private_exon_sequences.tsv" Outgroup Ceratotherium_simum

	sed -i 's/=OS//g' ${outdir}/*.fasta
	sed -i "s|GN=[^ ]*|GN=$goi|" ${outdir}/*.fasta

	done < "$gois"