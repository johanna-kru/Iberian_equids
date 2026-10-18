 aa_alignment <- read.alignment("/scratch_isilon/groups/compgen/jkrueger/horse-project/fastas/concatenation/concatenation_all_manual_trimmed_noX.fasta", seqtype="AA")

aa_dist <- dist.alignment(aa_alignment, matrix="similarity", gap=FALSE) # Fitch matrix for similarity

pca <- prcomp(as.matrix(aa_dist))

write.table(as.matrix(pca$x), file = "/scratch_isilon/groups/compgen/jkrueger/horse-project/fastas/concatenation/PCA_x_seqinr.txt", col.names=TRUE, quote=FALSE)


# Summary of the PCA
summary_pca <- summary(pca)

# Proportion of Variance
prop_variance <- summary_pca$sdev^2 / sum(summary_pca$sdev^2)

write.table(as.matrix(prop_variance), file = "/scratch_isilon/groups/compgen/jkrueger/horse-project/fastas/concatenation/prop_variance.txt", quote=FALSE)


