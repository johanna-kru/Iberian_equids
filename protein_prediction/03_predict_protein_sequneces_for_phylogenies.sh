#!/bin/bash


###################################################################################################################
###################################################################################################################

# Last modified by JOHANNA KRUEGER

### LOW QUALITY SITES IN VCF WILL BE MASKED AS "X" IN PROTEIN!

# This code masks any position indicated with a bed file in argument 7
# This code outputs the two alelles when the gene is heterozygous, and just one alelle when it is homozygous. In the 
# latter case, the option used calls all the alternative nucleotide when a snp is found.

# The code not yet makes any filtering of the vcf variants, therefore the vcf data used should be previously filtered
# in order to get biologically relevant data. Working on implementing the filters.

# Usage:

# ./thisprogram.sh gtf_file_location.gtf ref_genome_location.fa vcf_file.vcf sample_name output/abs/dir Species_name

# Three files are going to be written to the directory you give as the fifth argument in the "Usage". One two fasta files,
# one with the protein sequence of the input genes, and one with the correspondent DNA sequence. The third file is an intermediate
# tsv extension file with the CDS's retreived from the GTF file.

# IMPORTANT: The version of this program expects the name
# of the species whith the one you are working as a sixth argument in the initial command. This name will be appended
# to the fasta headers in the ouput.

# In order to translate proteins, the correspondent Combine_CDS_concatenator.py file must be in the folder where
# you are running this bash program.


# As stated in the usage, the fifth is the output directory. It is very important for it to be an absolute path, or a path
# that is contained in your current working directory.

###################################################################################################################
###################################################################################################################

#module load gcc/6.3.0
#module load xz/5.2.2
#module load BCFTOOLS/1.9
#module load SAMTOOLS/1.9
#module load zlib/1.2.4
#module load PYTHON/3.8.6

echo

gtf=$1
fasta=$2
vcf=$3
sample=$4
base_output=$5 
species=$6
bed=$7

gene_ids=$(cat $gtf | awk '{print $10}'  | sed -e 's/\;//g' -e 's/\"//g' | sort | uniq)


# Getting the sample name in useful form, specifically from Ludovic's horses sample naming
# If a simple sample name is provided, it is respected as-it-is
if [[ $sample == */* ]]; then 

	name=$(cut -d '/' -f 7 <<< $sample)


	id=$(cut -d '.' -f 1 <<< $name)"."$(cut -d '.' -f 2 <<< $name)

else id=$sample 

fi

rm -f ${base_output}/${id}.private_exon_sequences.tsv 2> /dev/null


for identifier in ${gene_ids[@]}; do

	# Checking if the gene in the loop is heterozygous or not

	Het_check=$( cat $gtf | grep $identifier | grep CDS | sort -nk4 | \
	while read l;do 
    
		scaffold=$(echo $l | awk '{print $1}');
		start_co=$(echo $l| awk '{print $4}');
		stop_co=$(echo $l | awk '{print $5}');

		genotype=$(bcftools view -H $vcf --samples $sample ${scaffold}:${start_co}-${stop_co} | \
		awk '{print $10}' | cut -d ':' -f 1)

		for i in $genotype; do if [[ $i == "0/1" ]]; then echo 1;
			 
			fi; 
		done;
	done; )

	
	if [ -z "$Het_check" ]; then

#############################################################################################################

		# This code is run if gene is homozygous

		echo $sample is homozygous for $identifier, getting private exons > /dev/stderr
  
		cat $gtf | grep $identifier | grep CDS | sort -nk4 | while read l;do 
    
			scaffold=$(echo $l | awk '{print $1}');
			start_co=$(echo $l| awk '{print $4}');
			stop_co=$(echo $l | awk '{print $5}');
			orientation=$(echo $l | awk '{print $7}');
			gene=$(echo $l | awk '{print $10}' | sed  -e 's/\;//g' -e 's/\"//g');
			tag=$(echo $l|awk '{print $12}' | sed -e 's/\;//g' -e 's/\"//g')"_A1";

			samtools faidx $fasta ${scaffold}:${start_co}-${stop_co} | \
			bcftools consensus --sample $sample -i "GT='AA'" $vcf --mask $bed | \
			bioawk -cfastx \
				-vtag=$tag \
				-vstop_co=$stop_co \
				-vstart_co=$start_co \
				-vscaffold=$scaffold\
				-vorientation=$orientation\
				-vgene=$gene\
			'{OFS="\t";print scaffold, start_co-1, stop_co, orientation, tag, $seq, gene}' | \
			sort -k1,1 -k2,2n >> ${base_output}/$id.private_exon_sequences.tsv;
		done;

	else

#############################################################################################################

		# This code is run if gene is heterozygous

		echo $sample is heterozygous for $identifier, getting private exons > /dev/stderr
  
		cat $gtf | grep $identifier | grep CDS | sort -nk4 | while read l;do 
    
			scaffold=$(echo $l | awk '{print $1}');
			start_co=$(echo $l| awk '{print $4}');
			stop_co=$(echo $l | awk '{print $5}');
			orientation=$(echo $l | awk '{print $7}');
			gene=$(echo $l | awk '{print $10}' | sed  -e 's/\;//g' -e 's/\"//g');
			tag=$(echo $l|awk '{print $12}' | sed -e 's/\;//g' -e 's/\"//g')"_A2";

			samtools faidx $fasta ${scaffold}:${start_co}-${stop_co} | \
			bcftools consensus --sample $sample -H 1 $vcf --mask $bed | \
			bioawk -cfastx \
				-vtag=$tag \
				-vstop_co=$stop_co \
				-vstart_co=$start_co \
				-vscaffold=$scaffold\
				-vorientation=$orientation\
				-vgene=$gene\
			'{OFS="\t";print scaffold, start_co-1, stop_co, orientation, tag, $seq, gene}' | \
			sort -k1,1 -k2,2n >> ${base_output}/$id.private_exon_sequences.tsv;
		done;
  
		cat $gtf | grep $identifier | grep CDS | sort -nk4 | while read l;do 
    
			scaffold=$(echo $l | awk '{print $1}');
			start_co=$(echo $l| awk '{print $4}');
			stop_co=$(echo $l | awk '{print $5}');
			orientation=$(echo $l | awk '{print $7}');
			gene=$(echo $l | awk '{print $10}' | sed  -e 's/\;//g' -e 's/\"//g');
			tag=$(echo $l|awk '{print $12}' | sed -e 's/\;//g' -e 's/\"//g')"_A1";


			samtools faidx $fasta ${scaffold}:${start_co}-${stop_co} | \
			bcftools consensus --sample $sample -H 2 $vcf --mask $bed | \
			bioawk -cfastx \
				-vtag=$tag \
				-vstop_co=$stop_co \
				-vstart_co=$start_co \
				-vscaffold=$scaffold\
				-vorientation=$orientation\
				-vgene=$gene\
			'{OFS="\t";print scaffold, start_co-1, stop_co, orientation, tag, $seq, gene}' | \
			sort -k1,1 -k2,2n >> ${base_output}/$id.private_exon_sequences.tsv;
		done;

#############################################################################################################

	fi;

done;


echo combining exons of provided gtf, genes translated are: $gene_ids > /dev/stderr


######
### Polishing output
######

###rename output because too many dots cause errors in parsind in next python script -> changed back because script did not work otherwise, not sure why
rename .realigned.q25.bam _realigned_q25_bam "${base_output}/$id.private_exon_sequences.tsv"


#### eliminate $id so that it parses renamed file
python3 Combine_CDS_concatenator.py "${base_output}/$id.private_exon_sequences.tsv" $id $species

