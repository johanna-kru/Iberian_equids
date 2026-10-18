 
  ### extract sequences matching headers
  fasta="$1"      #fasta
  headers="$2"    # headers without ">"

label=$(basename "$headers" .txt)
label=${label#headers_}


  awk -v var="$headers" 'BEGIN{while((getline<var)>0)l[">"$1]=1}/^>/{f=l[$1]}f' "$fasta" > "selected_"${label}".fasta"
