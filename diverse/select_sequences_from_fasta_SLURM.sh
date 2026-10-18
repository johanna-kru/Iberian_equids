#!/bin/bash

#SBATCH --time=00:10:00
#SBATCH --job-name=select_Seq
#SBATCH --mem-per-cpu 8000
#SBATCH -c 1
#SBATCH --partition compute



### created with help of a large language model

######
### Arguments
######

indir=$1 # gene names must be in uniprot format GN=GENE
outdir=$2
key=$3   # list with keyword for headers to keep (string or file)
filenames=$4
file=$(awk -v num=$SLURM_ARRAY_TASK_ID 'NR==num' $filenames)

######
### Prepare run
######

mkdir -p "$outdir"

######
### Extract sequences of interest from fasta
######

echo "Extracting sequences for $file..."

### Determine how to use the key: file or string
if [[ -f "$key" ]]; then
    grep_option=(-f "$key")
else
    grep_option=(-e "$key")
fi


### find input files
        # Enable nullglob so unmatched patterns disappear
    shopt -s nullglob

    # Expand the glob into an array
    db=( "$indir"/*"$file"*.fasta )

    # Check if any files were found
    if [ ${#db[@]} -eq 0 ]; then
        echo "No files found matching pattern: $indir/*$file*.fasta"
        exit 1
    fi

    # Print files for debugging
    echo "Files found: ${db[@]}"


    ### extract headers matching key(s)
    grep ">" "${db[@]}" | grep "${grep_option[@]}" | cut -d">" -f2 > "$outdir/headers_${file}.tmp"


    headers="$outdir/headers_${file}.tmp"

    ### extract sequences matching headers
  awk -v var="$outdir/headers_${file}.tmp" 'BEGIN{while((getline<var)>0)l[">"$1]=1}/^>/{f=l[$1]}f' "${db[@]}" > "$outdir/${file}_selected.fasta"
    rm "$headers"
