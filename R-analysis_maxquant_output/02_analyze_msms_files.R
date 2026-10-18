rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggridges)
library(RColorBrewer) # colors 
library("viridis") # colors 

######
### Read in data
######

setwd("Ranalysis_RUN9")

### list of filenames

fm <- read.table(file = "./msms_all/filenames.txt", stringsAsFactors = FALSE)$V1

### link between IDs and gene names
id2gn <- read.table(file = "id2genename_SAP_test_db_BCN_horses.txt", stringsAsFactors = FALSE)
colnames(id2gn) <- c( "protein", "id")
id2gn <- id2gn[!duplicated(id2gn), ]


ages <- read.csv(file = "samples_ages.csv", sep ="\t" , stringsAsFactors = FALSE)

### define samples
samples <-  sub(".txt", "", sub("msms_", "", fm) ) 

### master colors
mpal <- c("#2853a8ff","#7E97CBFF", 
          "#c32f31ff", "#FB8072" , 
          "#2f8450ff", "#62c1b8ff" , 
          "#7b43a9ff","#b98ed7ff", 
          "chocolate3","#FDB462" )


### colors for plot in supplement
pal <- c("#FDBF6F", "#E31A1C", "#FB9A99", "darkgreen", "#A6CEE3" , "#1F78B4",  "mediumpurple", "#FFED6F", "#CCEBC5", "#BC80BD", "#FCCDE5",
         "#B3DE69" , "#FDB462", "#80B1D3", "#FB8072", "#BEBADA", "#FFFFB3", "#8DD3C7", "#D9D9D9" )


### colors to highlight samples with long peptides
colsdens <-  c("#E31A1C", "#D9D9D9", "#D9D9D9", "#D9D9D9", "#D9D9D9", "#D9D9D9", "#D9D9D9",
               "#D9D9D9", "#D9D9D9", "#D9D9D9", "#D9D9D9", "#D9D9D9", 
               "#D9D9D9", "#338CC8", "#FDB462", "#69478bff", "#D9D9D9","#D9D9D9",
               "#D9D9D9", "#D9D9D9", "#D9D9D9", "#D9D9D9", "#D9D9D9" )


######
### Stats with accepted scans only
######

dir.create("figures")
dir.create("tables")



### RUN ic-score.R to get scores for each scan

source("ic-score.R")


######
### Get proteins with accepted peptides
######

prots <- list()

# filter by minimum ric score (relaxed score, 0.8 best score for shortest peptide of 7 aa without b1 nor y1)
minric=0.9

print( paste("Minimal ric score = ", minric, sep="") )

for( ii in fm ){
  
  pp <- msms_sel[[ii]][msms_sel[[ii]]$ricScore >= minric, "protein"]
  
  #  pp <- unique( sort( pp ) )
  pp <- sort( pp ) 
  prots[[ii]] <- pp
}

all_prots <- sort( unlist(prots) ) 

df_all_prots <- data.frame(
  protein = unique(all_prots),
  count = as.numeric(table(all_prots))
)

df_all_prots <- df_all_prots[!grepl("^CON_", df_all_prots[, 1]), ]
df_min_prots <- df_all_prots[df_all_prots$count != 1, ]

system('mkdir -p ./tables')

write.csv( df_min_prots, file="./tables/proteins_min2_accepted_peptides.csv", row.names = FALSE )




### Create table of proteins found
totalprot <- c()

samtot <- length(fm)


for(ii in fm){
  
  print(ii)
  sample <-  sub(".txt", "", sub("msms_", "", ii) ) 
  prots <- unique( sort( msms_sel[[ii]]$protein ) )
  
  totalprot <- c(totalprot, prots )
}

totalprot <- sort(totalprot)

prottab <- as.data.frame( table(totalprot) )
colnames(prottab) <- c("protein", "samples")

poitab <- prottab
poitab <- poitab[ !grepl("^CON_", poitab$protein ) , ]
poitab <- poitab[ !grepl("^REV_", poitab$protein ) , ]

poitab$nondis <- samtot-poitab$samples

poitab$protein <- factor(poitab$protein, levels = poitab$protein[order(poitab$samples)])
poitab <- melt(poitab)

### Evidence of accepted msms counts (df_min_prots from ic-score.R) by protein    
prot_evidence <- merge( x = poitab, y = df_min_prots, all = TRUE , by.x = "protein", by.y = "protein")
prot_evidence <- prot_evidence[ prot_evidence$variable == "samples" &  is.na(prot_evidence$variable) == FALSE , -2 ]
colnames(prot_evidence) <- c("protein", "samples_containing_protein_count", "accepted_msms_count" )
prot_evidence <- prot_evidence[ rev( order( prot_evidence$accepted_msms_count) ), ]
prot_evidence <- prot_evidence[ is.na( prot_evidence$accepted_msms_count) == FALSE, ]
dim(prot_evidence)
#View(prot_evidence)

write.csv(prot_evidence, file = "./tables/evidence_accepted_proteins_and_msms.csv", row.names = FALSE)

#######################################################################################
### For early iterations of MQ runs  
### isolate peptides of proteins that need manual check for low-complexity peptides

#  lc <- c("ZFC3H1", "MGP", "WAS", "PRELP", "XIRP2", "TTR", "TFIP11", "TGFB2", "DSPP", "KLK4", "TF", "NGEF", "DMP1", "F2", "OMD", "BGN", "BGLAP", "CLU", "SPARC", "KLKB1", "SERPINF1", "SPP1")

#  sub_msms <- list()
#  df_sub_msms <- data.frame()


# for(ii in fm){
#     
#     filn <- rep(ii , sum( msms_sel[[ii]]$protein %in% lc, na.rm = TRUE) ) # sum with na.rm counts only TRUE values in boolean vector
#     sub_msms[[ii]] <- msms_sel[[ii]][ msms_sel[[ii]]$protein %in% lc , ]
#     sub_msms[[ii]] <- cbind(sub_msms[[ii]], filn)
#     
#     df_sub_msms <- rbind(df_sub_msms, sub_msms[[ii]] )
#   
#   }
#  rm(sub_msms)
#  df_sub_msms <- df_sub_msms[order(df_sub_msms$protein), c("Sequence", "protein", "filn", "Scan.number", "PEP", "Score", "icScore", "ricScore")]
#  df_sub_msms <- df_sub_msms[df_sub_msms$ricScore >= 0.8, ]
#  dim(df_sub_msms)

#  write.csv(df_sub_msms, file = "./tables/scans_check_for_low_complexity.csv", row.names = FALSE )
####################################################################################  


### Create aa count tables only from accepted peptides
minric=1
msms_cutoff <- list()
msms_cutoff_all <- list()
msms_cutoff_count <- data.frame()
dir.create("fastas")

for( ii in fm ){
  
  print(ii)
  sample <-  sub(".txt", "", sub("msms_", "", ii) ) 
  
  msms_cutoff_all[[ii]] <- msms_sel[[ii]][ msms_sel[[ii]]$ricScore >= minric , ]
  msms_cutoff[[ii]] <- msms_sel[[ii]][ msms_sel[[ii]]$ricScore >= minric & !grepl('^CON_', msms_sel[[ii]]$mainID) & !grepl('^REV_', msms_sel[[ii]]$mainID), ]
  

  
  # count accepted spectra and add to data frame
  addf <- as.data.frame( table(msms_cutoff_all[[ii]]$protein) )
  try( addf$samples <- sample ) # use try to skip empty entries of msms_cutoff
  msms_cutoff_count <- rbind(msms_cutoff_count, addf)
  
  # create fasta file  
  fasta_file <- file( paste("./fastas/", sample, "_accepted_peptides.fasta", sep="") , "w")
  
  # Iterate over each row in the data frame
  for (jj in 1:nrow( msms_cutoff[[ii]] )) {
    # Construct the FASTA entry
    fasta_entry <- paste0(">", msms_cutoff_all[[ii]]$mainID[jj], " id=", msms_cutoff_all[[ii]]$id[jj], " PeptideID=", msms_cutoff_all[[ii]]$Peptide.ID[jj], " AndromedaScore=", msms_cutoff_all[[ii]]$Score[jj], " PEP=", msms_cutoff_all[[ii]]$PEP[jj], " bonds=" , msms_cutoff_all[[ii]]$bonds[jj] , " icScore=", msms_cutoff_all[[ii]]$icScore[jj], " ricScore=", msms_cutoff_all[[ii]]$ricScore[jj],  " bionScore=", msms_cutoff_all[[ii]]$bionScore[jj], " yionScore=", msms_cutoff_all[[ii]]$yionScore[jj], " GN=", msms_cutoff_all[[ii]]$protein[jj],         "\n", msms_cutoff_all[[ii]]$Sequence[jj], "\n")
    
    # Write the entry to the file
    cat(fasta_entry, file = fasta_file, append = TRUE)
  }
  
  # Close the file connection
  close(fasta_file)
  
}

######
### Get IDs of leading razor proteins. 
######

### Every protein of which a unique peptide (unique for the protein, i.e. this particular isoform) is found, will be the leading razor of its own protein group in MaxQuant  

# combine samples
totalmsms <- data.frame()

for(ii in fm){
  
  print(ii)
  sample <-  sub(".txt", "", sub("msms_", "", ii) ) 
  
  ### Eliminate rows with contaminants and reverse hits 
  wona <- !is.na(msms_cutoff[[ii]]$protein)
  addf <- msms_cutoff[[ii]][ wona , c("protein", "Sequence", "Length", "mainID") ] 
  try( addf$samples <- sample ) # use try to skip empty entries of msms_cutoff
  totalmsms <- rbind(totalmsms, addf)
  
}

# subset to samples produced in this study
firstsam <-  unique( sort(totalmsms[ grepl("P.*A", totalmsms$samples , perl = TRUE )  , "samples"] ) )

# exchange CPH sample (P1327A) for BCN version of same individual (P1327B)
firstsam[firstsam == "P1327A"] <- "P1327B"

# E.hydruntinus specimens (teeth) from same individual, compare only one sample from one tooth per individual
firstsam<- firstsam[ !( firstsam %in% c("P2856A", "P2857A") ) ] 

acceptedIDFirst <- totalmsms[totalmsms$samples %in% firstsam  , "mainID"]
acceptedIDFirst <- unique( sort(acceptedIDFirst)  )

wocon <- !grepl( paste( c("CON_", "REV_"), collapse="|"), acceptedIDFirst)
acceptedIDFirst <- acceptedIDFirst[ wocon]

write.table(acceptedIDFirst, file = "./tables/acceptedmsms_MainIds.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


######
### Plot peptide lengths
######

### 1. all by sample 

# combine samples
totalmsms_all <- data.frame()

for(ii in fm){
  
  print(ii)
  sample <-  sub(".txt", "", sub("msms_", "", ii) ) 
  
  ### Eliminate rows with contaminants and reverse hits 
  wona <- !is.na(msms_cutoff_all[[ii]]$protein)
  addf <- msms_cutoff_all[[ii]][ wona , c("protein", "Sequence", "Length", "mainID") ] 
  try( addf$samples <- sample ) # use try to skip empty entries of msms_cutoff
  totalmsms_all <- try(rbind(totalmsms_all, addf))
  
}

# only samples from this study
totalmsms_all <- totalmsms_all[totalmsms_all$samples %in% firstsam,]
totalmsms_all <- merge( totalmsms_all , ages, by="samples", all = FALSE)
totalmsms_all$samples <- factor(totalmsms_all$samples, levels = unique( totalmsms_all$samples[order(totalmsms_all$ageBP_Ma)]) )
# summarize contaminants
totalmsms_all[ grepl( paste( c("CON_"), collapse="|"), totalmsms_all$protein ), "protein"] <- "possible_contaminant"
# order levels of proteins with possible contaminant last
totalmsms_all$protein <- factor(totalmsms_all$protein, 
                                levels = c( unique(sort(as.character(totalmsms_all$protein)) )[unique(sort(as.character(totalmsms_all$protein) )) != "possible_contaminant"],
                                            "possible_contaminant") 
)



# density all in one plot

ggplot(totalmsms_all, aes(x = Length, color = samples)) +
  ylim(0, 0.15) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "lightgrey")) +
  geom_density(linewidth= 0.7, linetype = 1) +
  scale_color_manual(values = colsdens) + # Set color values
  geom_density(data = subset(totalmsms_all, samples == 'P2848A'), linewidth = 1.4)

width = 8
ggsave("./figures/Lengths_density_ACCEPTED_msms_all.pdf", width = width, height = 0.8 * width)


# density all ridgeline plot

# Calculate mean Length for each sample group
means_df <- totalmsms_all %>%
  group_by(samples) %>%
  summarise(mean_length = mean(Length)) %>%
  mutate(y = as.numeric(factor(samples, levels = rev(levels(totalmsms_all$samples)))))



ggplot(totalmsms_all, aes(x = Length, y = samples, fill = samples)) +
  geom_density_ridges(scale = 0.9, color=NA, fill = mpal[2], alpha = 1.0) +
  theme_ridges(line_size = 0.5, center_axis_labels = TRUE) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + 
  xlab("Length") +
  scale_x_continuous(limits = c(NA, 20)) +
  scale_y_discrete(limits = rev(levels(totalmsms_all$samples)))+  # needed to reverse order of factors
  geom_segment(data = means_df,
               aes(x = mean_length, xend = mean_length, y = y , yend = y + 0.6),
               color = mpal[3], size = 1)

width = 3
ggsave("./figures/Lengths_density_ACCEPTED_msms_ridges.pdf", width = width, height = 2 * width)




ggplot(totalmsms_all, aes(x = samples, y = Length, fill = samples)) +
  geom_violin(trim = FALSE, color = NA, fill = mpal[2], alpha = 1.0) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black"),
    legend.position = "none"
  ) +
  ylab("Length") +
  xlab(NULL) +
  scale_y_continuous(limits = c(NA, 20)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4,
               color = mpal[3])

width = 6
ggsave("./figures/Lengths_violin_ACCEPTED_msms.pdf",
       width = width, height = 0.5 * width)





# density by protein

ggplot(totalmsms_all, aes(x = Length, color = samples)) +
  ylim(0, 0.15) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "lightgrey"),
        panel.border = element_rect(colour = "lightgrey", fill=NA, linewidth=1.4) ) +
  facet_grid(protein ~ .) +
  geom_density(linewidth= 0.7, linetype = 1) +
  scale_color_manual(values = colsdens) + # Set color values
  geom_density(data = subset(totalmsms_all, samples == 'P2848A'), linewidth = 1.4)

width = 6
ggsave("./figures/Lengths_density_ACCEPTED_msms_by_protein.pdf", width = width, height = 2.8 * width)




######  
### Plot accepted scans
######

### Prepare data for plots
colnames(msms_cutoff_count) <- c("protein", "accepted_msms", "samples")

#proteins are noted as factors, so "possible_contaminant" must be added as possible factor
levels(msms_cutoff_count$protein) <- c(levels(msms_cutoff_count$protein), "possible_contaminant")

# summarize contaminants
msms_cutoff_count[ grepl( paste( c("CON_"), collapse="|"), msms_cutoff_count$protein ), "protein"] <- "possible_contaminant"

# aggregate by sample and by protein
msms_cutoff_count <- msms_cutoff_count %>%
  group_by(samples, protein) %>%
  summarize(accepted_msms = sum(accepted_msms))

# add ages
msms_cutoff_count <- merge(msms_cutoff_count, ages, by = "samples", all = TRUE)
msms_cutoff_count$samples <- factor(msms_cutoff_count$samples, levels = unique( msms_cutoff_count$samples[order(msms_cutoff_count$ageBP_Ma)] )  )

# order levels of proteins with possible contaminant last
msms_cutoff_count$protein <- factor(msms_cutoff_count$protein, 
                                    levels = c( unique(sort(as.character(msms_cutoff_count$protein)) )[unique(sort(as.character(msms_cutoff_count$protein) )) != "possible_contaminant"],
                                                "possible_contaminant") 
)



### count accepted msms

msms_cutoff_count_replicates <- msms_cutoff_count[grep("^P", msms_cutoff_count$samples), ]
write.csv(msms_cutoff_count_replicates, file="./tables/accepted_PSMs_replicates.csv")

msms_cutoff_count_First <- msms_cutoff_count[ msms_cutoff_count$samples %in% firstsam ,]
#msms_cutoff_count_First_sum <- summarise_all(group_by(msms_cutoff_count_First[, c(1,3)], samples), sum) # otherwise fill in ggplot will not work


ggplot(msms_cutoff_count_First, aes(x = samples, y = accepted_msms)) +
  geom_bar(stat='identity', position= position_dodge(0.7), width=0.7, fill = mpal[1]) + 
  ylab("accepted msms scans") +
  scale_x_discrete(limits = rev(levels(totalmsms_all$samples))) + # needed to reverse order of factors
  coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13)) +
  theme(axis.text.y = element_text(size = 13)) +
  theme(legend.position="none")


width = 12
ggsave("./figures/count_ACCEPTED_msms.pdf", width = width, height = 0.6 * width)



### relative frequencies accepted msms

ggplot(msms_cutoff_count_First, aes(x = samples, y = accepted_msms, fill = protein)) +
  geom_bar(stat='identity', position="fill", width=0.7, ) +
  scale_x_discrete(limits = rev(levels(totalmsms_all$samples))) + # needed to reverse order of factors
  coord_flip() +
  scale_fill_manual(values = rep(pal, length.out=19)) +
  scale_y_reverse() +
  ylab("fraction accepted msms scans") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13)) +
  theme(axis.text.y = element_text(size = 13)) 


width = 12
ggsave("./figures/rel_ACCEPTED_msms.pdf", width = width, height = 0.6 * width)


### Correlation sample diversity with sample age
shannon_by_sample <- msms_cutoff_count_First %>%
  group_by(samples) %>%
  summarise(
    ageBP_Ma = first(ageBP_Ma),
    total_msms = sum(accepted_msms, na.rm = TRUE),
    n_proteins = sum(accepted_msms > 0, na.rm = TRUE),
    # pointed brackets indicate that only last expression is returned
    shannon_H = {
      counts <- accepted_msms[!is.na(accepted_msms) & accepted_msms > 0]
      p <- counts / sum(counts)
      -sum(p * log(p))
    },
    .groups = "drop" # grouping removed, may be useful for subsequent calculations on output matrix
  )


# Pearson correlation
cor.test(
  shannon_by_sample$ageBP_Ma,
  shannon_by_sample$shannon_H,
  method = "pearson"
)

# data:  shannon_by_sample$ageBP_Ma and shannon_by_sample$shannon_H
# t = -6.9041, df = 18, p-value = 1.868e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.9400798 -0.6572820
# sample estimates:
#   cor 
# -0.8519917 

# Spearman correlation
cor.test(
  shannon_by_sample$ageBP_Ma,
  shannon_by_sample$shannon_H,
  method = "spearman",
  exact = FALSE
)

# data:  shannon_by_sample$ageBP_Ma and shannon_by_sample$shannon_H
# S = 2028.7, p-value = 0.01738
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#   rho 
# -0.5253286 
