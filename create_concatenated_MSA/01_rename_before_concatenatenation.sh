#!/bin/bash


indir="ASRs_to_ref_manual"

outdir="ASRs_to_ref_manual/renamed"

gois="ASRs_to_ref_manual/gois_to_rename_and_concatenate.txt" #ls *.fasta | cut -d"_" -f1 | grep -v -E "AMELY|COL2A1|SERPINA1" > gois_to_rename_and_concatenate.txt


### Create outdir 

mkdir -p $outdir


### Rename headers, so that sequences can be concatenated

while read goi; do

    cat $indir/*$goi*".fasta" | sed 's/>|Outgroup/>white_rhino|Outgroup/g' | awk -F'_' '/consensus/ {print ">"$2"_consensus|"; next} {print}' > $outdir/$goi"_renamed.prefasta1"
    # add name for rhino headers before | (cutting delimiter), remove gene name before utting delimiter for ASRs

    cat $outdir/$goi"_renamed.prefasta1" | awk -F "|" '/^>/ {print $1; next} 1' | sed 's/\r//' > $outdir/$goi"_renamed.prefasta2"
    # sed 's/\r//' converts line endings from windows to linux format


    awk '/^>/{if($0 in seen){$0 = $0 "_" ++seen[$0]} else {seen[$0]}}1' $outdir/$goi"_renamed.prefasta2" > $outdir/$goi"_renamed.fasta"
    # appends underscore and running numbering to duplicate headers


    rm $outdir/$goi"_renamed.prefasta"*

done < "$gois"


