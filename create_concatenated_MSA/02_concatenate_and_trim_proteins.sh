#!/bin/bash
### ran as interactive job

### all individuals should have same fasta header across all input MSAs!

indir="renamed"
outdir=$indir

######
### Concatenate and remove X (optional)
######


catfasta2phyml.pl --fasta --concatenate --verbose $indir/*renamed.fasta \
> $outdir/concatenated_aln.prefasta1


######
### Trim lines
######

/bin/trimAl/source/trimal -in concatenated_aln.prefasta1 -out $outdir/concatenated_aln.fasta -htmlout $outdir/concatenated_aln.fasta.html -resoverlap 0.75 -seqoverlap 80

