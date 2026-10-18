#!/bin/bash


### Checks database for completeness and errors in headers

######
### Arguments
####


db=$1




######
### Get stats
#####


echo "$db"

echo 'Check number of ">": cat $db | grep ">" | wc -l'
cat $db | grep ">" | wc -l

echo 'Check number of headers: cat $db | grep "^>" | wc -l'
cat $db | grep "^>" | wc -l

echo 'Check if every entry has "GN=": cat $db | grep "^>" |  grep "GN=" | wc -l'
cat $db | grep "^>" | grep "GN=" | wc -l

echo 'Check if every entry has "OS=": cat $db | grep "^>" | grep "OS=" | wc -l'
cat $db | grep "^>" | grep "OS=" | wc -l

echo 'Check if every entry has "=SO": cat $db | grep "^>" | grep "=SO" | wc -l'
cat $db | grep "^>" | grep "=SO" | wc -l

echo 'Check if every entry has two pipes: cat $db | grep "^>" | grep -o "|" | wc -l'
cat $db | grep "^>" | grep -o "|" | wc -l

echo 'Check if database has too long identical beginning of headers, may cause problems in string matchig in awk: cat $db |  grep "^>" | grep -o "UniProtKB/Swiss-Prot proteins" | wc -l'
cat $db |  grep "^>" | grep -o "UniProtKB/Swiss-Prot proteins" | wc -l

echo 'Check for word "isoform" (may cause problems in MQ and GP): cat $db | grep "^>" | grep -o -i "isoform" | wc -l'
cat $db | grep "^>" | grep -o -i "isoform" | wc -l

echo 'Check for commas, which may disrupt output csvs: cat $db | grep "^>" | grep -o "," | wc -l'
cat $db | grep "^>" | grep -o "," | wc -l

echo 'Check for semicolons, which may disrupt output csvs: cat $db | grep "^>" | grep -o ";" | wc -l'
cat $db | grep "^>" | grep -o ";" | wc -l

echo 'Check for colons, which may disrupt output csvs: cat $db | grep "^>" | grep -o ":" | wc -l'
cat $db | grep "^>" | grep -o ":" | wc -l

echo 'Check for species in db cat $db | grep OS= | awk -F "OS=" {print $2} | awk -F" " {print $1, $2} | sort | uniq -c'
cat $db | grep OS= | awk -F "OS=" '{print $2}' | awk -F" " '{print $1, $2}' | sort | uniq -c

echo 'List of genes'
cat $db | grep ^">" | awk -F "GN=" '{ print $2}' | cut -d" " -f1 | sort | uniq -c | sort -nr

echo 'number of genes'
cat $db | grep ^">" | awk -F "GN=" '{ print $2}' | cut -d" " -f1 | sort | uniq -c | wc -l




