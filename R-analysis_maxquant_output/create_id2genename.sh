#!/bin/bash
#Create file that links gene names to uniprot and genebank IDs

db=$1

cat $db | grep ">" | awk -F "|" '{print $3,$2}' | awk -F 'OS=|=SO' '{print $0,$2}' | rev | cut -d" " -f1,2,3 | rev > id_to_genename.txt
cat $db | grep ">" | awk -F 'GN=' '{print $2}' | cut -d" " -f1 > genenames.txt

cat genenames.txt | wc -l


cat id_to_genename.txt | wc -l


#print files next to each other
#https://stackoverflow.com/questions/13341832/display-two-files-side-by-side
pr  -J -m -t genenames.txt id_to_genename.txt | awk '{print $0,$2,$3,$4}' |  awk '{print $1,$7}' > id2genename.txt


cat id2genename.txt | wc -l

rm id_to_genename.txt
rm genenames.txt
