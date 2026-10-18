
inaln="$1"

filestem=$(echo $inaln | awk '{sub(/\.fasta$/, ""); print}')
cat ${inaln} | sed '/^>/! s/X/-/g' > ${filestem}"_noX.fasta"