#!/bin/bash

# written by Johanna Krueger
# version 2024-04-11

# Script to algin databases by gene of interest, i.e. protein
# in every header the gene must be indicated by GN=<NAME>.

# Entries that do not fulful this condition will be ignored.

# If the headers are identical for the first 21-30 characters (I am not sure yet for how long they can be identical without problem),
# make sure to shorten the identical expression, otherwise awk will regard them as identical.

# If script runs withour problems, uncomment deleting of temporaty files at end of script.

# DEPENDENCIES: mafft (any version should so)
# USAGE: bash get_aln_from_db.sh <datbase.fasta> 
# OPTIONAL ARGUMENTS: bash get_aln_from_db.sh <datbase.fasta> <output/directory/here> <prefix_for_output_files> <list_of_genes.txt> <list_of_species.txt>


#####################################################################



######
### Arguments
######

#database
dbin="$1" # db: gene names must be in uniprot format GN=GENE

# [OPTIONAL] define output directory
outdir="$2"
if [ "$2" != "" ]; then
	outdir="$2"
else
	outdir="./"
fi
# [OPTIONAL] define string that is in all output file names

if [ "$3" != "" ]; then
	outfilestem="$3"
else
	outfilestem="aligned"

fi

# [OPTIONAL] txt with genes of interest that should be aligned, one gene per line
gois="$4"

# [OPTIONAL] species of interest in txt file, one species per line
if [ "$5" != "" ]; then
	spois="$5" 
else
	spois=""
fi


######
### Change database headers because some formats may cause downstream problems
######

### Shorten "UniProtKB/Swiss-Prot proteins", because otherwise awk will match string because the beginning of very long headers is the same

# If you get other proteins in your alignment, most likely their headers are similar in the 
# beginning to some headers of your protein that you want to align.
# Just shorten them similar to this example.

sed 's/UniProtKB\/Swiss-Prot proteins/UniProtKB\/Swiss-Prot/g' $dbin > $outdir/"tmp_db.fasta"

db=$outdir/"tmp_db.fasta"


######
### Get Gene of interest
######

if [ "$4" != "" ]; then
	gois="$4" 
else

  cat $db | grep ">" | grep "GN=" | awk -F "GN=" '{print $2}' | cut -d" " -f1 | sort | uniq > $outdir/"gois_generated_from_db.txt"

	gois=$outdir/"gois_generated_from_db.txt"
fi



while read goi; do


######
#### Prepare directories
#######

  mkdir -p $outdir/$goi"_tmp"


######
### Extract gois from database
######

  cat $db | grep ">" | grep "GN=$goi"  | cut -d">" -f2  > $outdir/$goi"_tmp/headers_${goi}.txt"
  
  headers=$outdir/$goi"_tmp/headers_${goi}.txt"
 
  awk -v var="$headers" 'BEGIN{while((getline<var)>0)l[">"$1]=1}/^>/{f=l[$1]}f' $db > $outdir/$goi"_tmp/all_${goi}_unaligned.fasta"

 
 
 # alternative with seqtk: seqtk subseq $db $outdir/$goi"_tmp"/headers_goi.txt > $outdir/$goi"_tmp"/all_goi_unaligned.fasta

######
### Subset to species of interest (OPTIONAL, only runs, if $5 is given)
#####

if [ "$5" != "" ]; then
  cat $outdir/$goi"_tmp/all_${goi}_unaligned.fasta" | grep ">" | grep -e -i -f "$spois"  | cut -d">" -f2  >  $outdir/$goi"_tmp"/headers_spoi.txt

  specs=$outdir/$goi"_tmp"/headers_spoi.txt

  awk -v var="$specs" 'BEGIN{while((getline<var)>0)l[">"$1]=1}/^>/{f=l[$1]}f' $outdir/$goi"_tmp/all_${goi}_unaligned.fasta" > $outdir/$goi"_tmp/${goi}_unaligned.fasta"
  
  # alternative with seqtk: seqtk subseq $outdir/$goi"_tmp"/all_goi_unaligned.fasta $outdir/$goi"_tmp"/headers_spoi.txt > $outdir/$goi"_tmp"/goi_unaligned.fasta

else
	mv $outdir/$goi"_tmp/all_${goi}_unaligned.fasta" $outdir/$goi"_tmp/${goi}_unaligned.fasta"

fi

######
### Align gois
######



  /scratch_isilon/groups/compgen/jkrueger/bin/mafft-7.490-without-extensions/bin/mafft \
    --6merpair \
    --maxiterate 1000 \
    --amino \
    --reorder \
    --bl 80 \
    $outdir/$goi"_tmp/${goi}_unaligned.fasta" > "$outdir/${outfilestem}_${goi}.fasta"
  # maxiterate: number cycles of iterative refinement are performed. Default: 0
  # bl: BLOSUM Matix, 80 because of conserved sequences


  ######
  ### Remove temporary files
  ######

 # rm -r $outdir/$goi"_tmp"

	done < "$gois"


 # rm $outdir/"tmp_db.fasta"
