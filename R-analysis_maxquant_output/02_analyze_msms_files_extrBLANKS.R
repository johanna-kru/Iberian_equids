rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggridges)
library(RColorBrewer) # colors 
library("viridis") # colors 

setwd("Ranalysis_extraction_blanks")

######
### Read in data
######

### list of filenames

fm <- read.table(file = "./msms_all/filenames.txt", stringsAsFactors = FALSE)$V1

### link between IDs and gene names
id2gn <- read.table(file = "id2genename_SAP_test_db_BCN_horses.txt", stringsAsFactors = FALSE)
colnames(id2gn) <- c( "protein", "id")
id2gn <- id2gn[!duplicated(id2gn), ]


ages <- read.csv(file = "samples_ages.csv", sep ="\t" , stringsAsFactors = FALSE)

### define samples
samples <-  sub(".txt", "", sub("msms_", "", fm) ) 

### get connection blank sample and biological samples
blk2sam <- read.table(file = "./evidence_all/associated_blanks.txt", stringsAsFactors = FALSE)
colnames(blk2sam) <- c("samples", "blank_name")

### define colors for most plots
pal <- c("#FDBF6F", "#E31A1C", "#FB9A99", "darkgreen", "#A6CEE3" , "#1F78B4",  "mediumpurple", "#FFED6F", "#CCEBC5", "#BC80BD", "#FCCDE5",
         "#B3DE69" , "#FDB462", "#80B1D3", "#FB8072", "#BEBADA", "#FFFFB3", "#8DD3C7", "#D9D9D9" )




######
### Stats with accepted scans only
######

dir.create("figures")
dir.create("tables")



### RUN ic-score.R to get scores for each scan

source("ic-score.R")

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
# In this case keep contaminants and reverse hits
#poitab <- poitab[ !grepl("^CON_", poitab$protein ) , ]
#poitab <- poitab[ !grepl("^REV_", poitab$protein ) , ]

poitab$nondis <- samtot-poitab$samples

poitab$protein <- factor(poitab$protein, levels = poitab$protein[order(poitab$samples)])
poitab <- melt(poitab)


#######################################################################################

### Create aa count tables only from accepted peptides
minric=0.8
msms_cutoff <- list()
msms_cutoff_count <- data.frame()
dir.create("fastas")

for( ii in fm ){
  
  print(ii)
  sample <-  sub(".txt", "", sub("msms_", "", ii) ) 
  
  msms_cutoff[[ii]] <- msms_sel[[ii]][msms_sel[[ii]]$ricScore >= minric, ]
  

  
  # count accepted spectra and add to data frame
  addf <- as.data.frame( table(msms_cutoff[[ii]]$protein) )
  try( addf$samples <- sample ) # use try to skip empty entries of msms_cutoff
  msms_cutoff_count <- rbind(msms_cutoff_count, addf)
  
  
  # create fasta file  
  fasta_file <- file( paste("./fastas/", sample, "_accepted_peptides.fasta", sep="") , "w")
  
  # Iterate over each row in the data frame
  for (jj in 1:nrow( msms_cutoff[[ii]] )) {
    # Construct the FASTA entry
    fasta_entry <- paste0(">", msms_cutoff[[ii]]$mainID[jj], " id=", msms_cutoff[[ii]]$id[jj], " PeptideID=", msms_cutoff[[ii]]$Peptide.ID[jj], " AndromedaScore=", msms_cutoff[[ii]]$Score[jj], " PEP=", msms_cutoff[[ii]]$PEP[jj], " bonds=" , msms_cutoff[[ii]]$bonds[jj] , " icScore=", msms_cutoff[[ii]]$icScore[jj], " ricScore=", msms_cutoff[[ii]]$ricScore[jj],  " bionScore=", msms_cutoff[[ii]]$bionScore[jj], " yionScore=", msms_cutoff[[ii]]$yionScore[jj], " GN=", msms_cutoff[[ii]]$protein[jj],         "\n", msms_cutoff[[ii]]$Sequence[jj], "\n")
    
    # Write the entry to the file
    cat(fasta_entry, file = fasta_file, append = TRUE)
  }
  
  # Close the file connection
  close(fasta_file)
  
}

######
### Get IDs of leading razor proteins. 
######

###Every protein of which a unique peptide (unique for the protein, i.e. this particular isoform) is found, qill be the leading razor of its own protein group in MaxQuant  

# combine information of msms:cutoff into one big data frame with all samples
totalmsms <- data.frame()

for(ii in fm){
  
  print(ii)
  sample <-  sub(".txt", "", sub("msms_", "", ii) ) 
  
  ### Eliminate rows with contaminants and reverse hits 
  wona <- !is.na(msms_cutoff[[ii]]$protein)
  addf <- msms_cutoff[[ii]][ wona , c("protein", "Sequence", "Length", "mainID") ] 
  try( addf$samples <- sample ) # use try to skip empty entries of msms_cutoff
  totalmsms <- try(rbind(totalmsms, addf))
  
}

# subset to samples produced in this study
firstsam <-  unique( sort(totalmsms[ grepl("20", totalmsms$samples , perl = TRUE )  , "samples"] ) )



acceptedIDFirst <- totalmsms[totalmsms$samples %in% firstsam  , "mainID"]
acceptedIDFirst <- unique( sort(acceptedIDFirst)  )

#wocon <- !grepl( paste( c("CON_", "REV_"), collapse="|"), acceptedIDFirst)
#acceptedIDFirst <- acceptedIDFirst[ wocon]

write.table(acceptedIDFirst, file = "./tables/acceptedmsms_MainIds.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


######
### Plot peptide lengths
######

### 1. all by sample 


totalmsms_all <- totalmsms

# only samples from this study
totalmsms_all <- totalmsms_all[totalmsms_all$samples %in% firstsam,]
totalmsms_all <- merge( totalmsms_all , ages, by="samples", all.x = TRUE)
totalmsms_all$samples <- factor(totalmsms_all$samples, levels = unique( totalmsms_all$samples[order(totalmsms_all$ageBP_Ma)]) )
totalmsms_all$protein <- factor(totalmsms_all$protein, 
                                levels = c( unique(sort(as.character(totalmsms_all$protein)) )[unique(sort(as.character(totalmsms_all$protein) )) != "possible_contaminant"],
                                            "possible_contaminant") 
)

totalmsms_all <- merge( totalmsms_all , blk2sam, by="samples", all.x = TRUE)
totalmsms_all$blank_name <- factor(totalmsms_all$blank_name)



# density all ridgeline plot

ggplot(totalmsms_all, aes(x = Length, y = blank_name, fill = blank_name)) +
  geom_density_ridges(scale = 1.9, fill = "#1F78B4", alpha = 0.5, color = "#1F78B4") +
  theme_ridges(line_size = 0.5, center_axis_labels = TRUE) +
  xlab("Length") +
  scale_y_discrete(limits = rev(levels(totalmsms_all$blank_name)))  # needed to reverse order of factors


width = 7
ggsave("./figures/Lengths_density_ACCEPTED_msms_ridges.pdf", width = width, height = 0.7 * width)



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
#msms_cutoff_count <- merge(msms_cutoff_count, ages, by = "samples", all = TRUE)
#msms_cutoff_count$samples <- factor(msms_cutoff_count$samples, levels = unique( msms_cutoff_count$samples[order(msms_cutoff_count$ageBP_Ma)] )  )

# add blank names
msms_cutoff_count <- merge(msms_cutoff_count, blk2sam, by = "samples", all.x = TRUE)
msms_cutoff_count$blank_name <- factor(msms_cutoff_count$blank_name)

# order levels of proteins with possible contaminant last
msms_cutoff_count$protein <- factor(msms_cutoff_count$protein, 
                                    levels = c( unique(sort(as.character(msms_cutoff_count$protein)) )[unique(sort(as.character(msms_cutoff_count$protein) )) != "possible_contaminant"],
                                                "possible_contaminant") 
)


### count accepted msms

msms_cutoff_count_First <- msms_cutoff_count[ msms_cutoff_count$samples %in% firstsam ,]


ggplot(msms_cutoff_count_First, aes(x = blank_name, y = accepted_msms, fill = blank_name)) +
  geom_bar(stat='identity', position= position_dodge(0.7), width=0.7, color = "#1b6aa2", fill= "#1b6aa2") + 
#  scale_fill_manual(values = bluehue) +
  ylab("accepted msms scans") +
  scale_x_discrete(limits = rev(levels(totalmsms_all$blank_name))) + # needed to reverse order of factors
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
# use ggplot2 to get relative frequencies right away without changing input table


ggplot(msms_cutoff_count_First, aes(x = blank_name, y = accepted_msms, fill = protein)) +
  geom_bar(stat='identity', position="fill", width=0.7, ) +
  scale_x_discrete(limits = rev(levels(totalmsms_all$blank_name))) + # needed to reverse order of factors
  coord_flip() +
  scale_fill_manual(values = rep(pal, length.out=19)) +
#  scale_y_reverse() +
  ylab("fraction accepted msms scans") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13)) +
  theme(axis.text.y = element_text(size = 13)) 


width = 12
ggsave("./figures/rel_ACCEPTED_msms.pdf", width = width, height = 0.6 * width)

