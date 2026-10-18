rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer) # colors 
library("viridis") # colors 

######
### Read in data
######

### list of filenames
fl <- read.table(file = "./proteinGroups_all/filenames.txt", stringsAsFactors = FALSE)$V1

### link between IDs and gene names
id2gn <- read.table(file = "id2genename_SAP_test_db_BCN_horses.txt", stringsAsFactors = FALSE)
colnames(id2gn) <- c( "protein", "id")
id2gn <- id2gn[!duplicated(id2gn), ]


ages <- read.csv(file = "samples_ages.csv", sep ="\t" , stringsAsFactors = FALSE)

### define samples
samples <-  sub(".txt", "", sub("proteinGroups_", "", fl) ) 

### define colors for most plots
pal <- c("#FDBF6F", "#E31A1C", "#FB9A99", "darkgreen", "#A6CEE3" , "#1F78B4",  "mediumpurple", "#FFED6F", "#CCEBC5", "#BC80BD", "#FCCDE5",
         "#B3DE69" , "#FDB462", "#80B1D3", "#FB8072", "#BEBADA", "#FFFFB3", "#8DD3C7", "#D9D9D9" )



### read in proteinGroups files
pgr <- list()
pgr_sel <- list()

for(ii in fl){
  
  print(ii)
  sample <-  sub(".txt", "", sub("proteinGroups_", "", ii) ) 
  
  # define relevant columns 
  relsel <- c( "mainID", paste( "Razor...unique.peptides.", sample, sep = "" ),  paste( "Unique.peptides.", sample, sep = "" ), "Unique...razor.sequence.coverage....", "Sequence.length" , "MS.MS.IDs")
  
  # use try so that loop continues if one file empty
  pgr[[ii]] <- try (read.csv( paste( "./proteinGroups_all/", ii, sep = "")  , sep = '\t', stringsAsFactors = FALSE) )
  
  # use grepl to get all columns containing "Protein.IDs" as the column name might differ between versions of MQ
  pgr[[ii]]$mainID <- pgr[[ii]][, grepl("Protein.IDs", colnames( pgr[[ii]] ))]
  
  # sapply function to predefined column of df
  pgr[[ii]]$mainID <- sapply( pgr[[ii]]$mainID, function(x)unlist(strsplit( as.character(x), ";" )) [1] )
  
  # subset to important columns
  pgr_sel[[ii]] <- pgr[[ii]][, names(pgr[[ii]]) %in% as.character( relsel ) ]
  pgr_sel[[ii]]$No.Scans <- sapply( pgr_sel[[ii]]$MS.MS.IDs , function(x)length( strsplit(as.character(x), split = ";")[[1]] ) )
  
  pgr_sel[[ii]][pgr_sel[[ii]] == 0] <- NA
  
  # Link IDs to protein names (without species, if species wanted add columns 3,4 in id2gn). 
  # REV and CON entries won't have a "protein" entry, so that "protein.long" will start with NA. This is needed to eliminate them later!
  pgr_sel[[ii]] <- merge(x = pgr_sel[[ii]], y = id2gn[, c(1,2)], all.x = TRUE, by.x = "mainID", by.y = "id")
  pgr_sel[[ii]]$protein.long <- paste( pgr_sel[[ii]]$protein, pgr_sel[[ii]]$mainID, sep = "-" ) 
  
  # Replace missing gene names with Protein.IDs (for contaminants)
  pgr_sel[[ii]]$protein[ is.na( pgr_sel[[ii]]$protein ) ] <- pgr_sel[[ii]]$mainID[ is.na( pgr_sel[[ii]]$protein ) ]
  
  # Calculate normalized coverage
  pgr_sel[[ii]]$Unique_razor_coverage_per_100aa <- round( pgr_sel[[ii]]$Unique...razor.sequence.coverage.... * 100 /  pgr_sel[[ii]]$Sequence.length , digits = 2 )
  
  # Calculate number of amino acids recovered
  pgr_sel[[ii]][, paste("aa_", sample, sep = "")] <- round( pgr_sel[[ii]]$Unique...razor.sequence.coverage.... * pgr_sel[[ii]]$Sequence.length / 100 , digits = 0 )
  
}

#View(pgr_sel$proteinGroups_P1327B.txt)
rm(pgr)






######
#### Plot characteristics per sample
######A

dir.create("figures")

cols <-c("indianred3", "deepskyblue4", "seagreen4", "lightsalmon3", "slateblue3"  )

### Function to plot protein groups against chosen characteristic (defined in type and colsel)

countplot <- function(intab, 
                      fn, 
                      id,
                      type = "razor + unique peptides",
                      rowsel = intab$protein.long[!grepl("^NA-", intab$protein.long )],  # eliminates contaminants and rev
                      colsel = c( "protein.long", paste( "Razor...unique.peptides.", id, sep = "") ), 
                      col = cols[1]) {
  
  intab <- intab[ intab$protein.long %in% rowsel ,  colnames(intab) %in%  colsel ]
  intab <-merge(intab, as.data.frame(rowsel), by.x = "protein.long", by.y = "rowsel", all = TRUE)
  # shorten long protein id
  intab[intab == "COL17A1-TP4_CGG_1_018394ilkProteins.MilkProteins_COL17A1.XM_014733101.2_A1"] <- "COL17A1-TP4_CGG_1_018394"
  mm <- melt(intab)
  colnames(mm) <- c("protein.long", type, "value")
  
  
  ggplot(mm, aes(x = protein.long, y = value)) +
    geom_bar(stat='identity', position= position_dodge(0.7), width=0.7, fill = col) +
    ylab(type) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
    theme(axis.text.x = element_text( angle = 45, hjust = 1) ) +
    theme(legend.position="none")
  
  
  width = 7
  ggsave(fn, width = width, height = 1 * width)
  
  return(intab)
}



######
### Plots all samples combined
######


#### 1. Plot total amount

totalaa <- as.data.frame(samples)
totalaa$total_aa <- NA

# sum up aa
for(ii in fl){
  
  print(ii)
  sample <-  sub(".txt", "", sub("proteinGroups_", "", ii) ) 
  
  # Eliminate rows with contaminats and reverse hits 
  wocon <- !grepl( paste( c("CON_", "REV_"), collapse="|"), pgr_sel[[ii]]$mainID )
  aa_wocon <- sum(pgr_sel[[ii]][ wocon , paste("aa_", sample, sep = "")], na.rm = TRUE)
  totalaa[totalaa$samples == sample , "total_aa"] <- aa_wocon
  
}




totalaa <- merge( totalaa, ages, by="samples", all = FALSE)
totalaa$samples <- factor(totalaa$samples, levels = unique( totalaa$samples[order(totalaa$ageBP_Ma)]) )

# subset to samples produced in this study
firstsam <-  totalaa[ grepl("P.*A", totalaa$samples , perl = TRUE )  , "samples"]

# exchange CPH sample (P1327A) for BCN version of same individual (P1327B)
firstsam[firstsam == "P1327A"] <- "P1327B"
firstsam<- firstsam[ !( firstsam %in% c("P2856A", "P2857A") ) ] # E.hydruntinus specimens (teeth) from same individual, compare only one sample from one tooth per individual
totalaaFirst <- totalaa[ totalaa$samples %in% firstsam  ,]

# add level by ages to samples so that sorted in plot
totalaaFirst$samples <- factor(totalaaFirst$samples, levels = unique( totalaaFirst$samples[order(totalaaFirst$ageBP_Ma)]) )



ggplot(totalaaFirst, aes(x = samples, y = total_aa)) +
  geom_bar(stat='identity', position= position_dodge(0.7), width=0.7, fill = "darkgreen") +
  ylab("amino acids recovered") +
  scale_x_discrete(limits = rev(levels(totalaaFirst$samples))) + # needed to reverse order of factors
  coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
  theme(axis.title = element_text(size = 28)) +
  theme(axis.text.x = element_text(size = 26)) +
  theme(axis.text.y = element_text(size = 26)) +
  theme(legend.position="none")


width = 12
ggsave("./figures/total_aa_protein_Groups.pdf", width = width, height = 1.6 * width)




### 2. keep aa by protein

totalprotls <- list()

for(ii in fl){
  
  print(ii)
  sample <-  sub(".txt", "", sub("proteinGroups_", "", ii) ) 
  
  # Eliminate rows with contaminants and reverse hits because also not counted in totalaa
  wocon <- !grepl( paste( c("CON_", "REV_"), collapse="|"), pgr_sel[[ii]]$mainID )
  prot_wocon <- pgr_sel[[ii]][ wocon , c( "protein.long",  paste("aa_", sample, sep = "") )  ]
  
  # otherwise empty dfs will break loops
  if (nrow(prot_wocon) != 0) {
    prot_wocon$samples <- sample
  }
  
  
  colnames(prot_wocon)[2] <- "aa"
  totalprotls[[ii]] <- prot_wocon
  
  
}

totalprot <- do.call(rbind, totalprotls)
totalprot$protein <- sapply(strsplit(totalprot$protein.long, "-"), "[[", 1)

# for each sample, group by protein and sum, up aa for each protein, alternative to aggregate
sum_totalprot <- totalprot %>%
  group_by(samples, protein) %>%
  summarize(total_aa = sum(aa))

sum_totalprot <- merge( sum_totalprot, ages, by="samples", all = FALSE)

# LEVELS MUST BE UNIQUE!
sum_totalprot$samples <- factor(sum_totalprot$samples, levels = unique( sum_totalprot$samples[order(sum_totalprot$ageBP_Ma)]) )

# subset to samples produced in this study
totalprotFirst <- sum_totalprot[  sum_totalprot$samples %in% firstsam ,]

# relative frequencies
relprotFirst <- merge(x = totalprotFirst, y = totalaaFirst, by="samples", all.x = TRUE)
relprotFirst$percprot <- relprotFirst$total_aa.x/relprotFirst$total_aa.y

ggplot(relprotFirst, aes(x = samples, y = percprot, fill = protein)) +
  geom_bar(stat='identity', position="stack", width=0.7, ) +
  scale_fill_manual(values = rep(pal, length.out=18)) +
  ylab("fraction amino acids recovered") +
  scale_x_discrete(limits = rev(levels(totalaaFirst$samples))) + # needed to reverse order of factors
  coord_flip() + #swap x and y
  scale_y_reverse() + # to make top entries from legend appear leftmost
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13)) +
  theme(axis.text.y = element_text(size = 13))

width = 12
ggsave("./figures/rel_prot_proteinGroups.pdf", width = width, height = 0.6 * width)


##### 2. Proteins recovered

### 2.1 proteins
totalprot <- c()

samtot <- length(fl)


for(ii in fl){
  
  print(ii)
  sample <-  sub(".txt", "", sub("proteinGroups_", "", ii) ) 
  prots <- unique( sort( pgr_sel[[ii]]$protein ) )
  
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


ggplot(poitab, aes(x = protein, y = value, fill=variable)) +
  geom_bar(stat='identity', position = position_fill(reverse = TRUE) ) +
  coord_flip() +
  ylab("fraction of samples") +
  scale_fill_manual("", values = c("darkgreen", "grey"), labels = c("detected", "not detected")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15)) +
  theme(axis.text.y = element_text(size = 15)) +
  theme(legend.position="bottom")


width = 12
ggsave("./figures/proteins_recovered_proteinGroups.pdf", width = width, height = 1.5 * width)


### 2.1 Isoforms recovered

totalprot <- c()


for(ii in fl){
  
  print(ii)
  sample <-  sub(".txt", "", sub("proteinGroups_", "", ii) ) 
  prots <- unique( sort( pgr_sel[[ii]]$protein.long ) )
  
  totalprot <- c(totalprot, prots )
}

totalprot <- sort(totalprot)

isotab <- as.data.frame( table(totalprot) )
colnames(isotab) <- c("protein", "samples")


isotab <- isotab[ !grepl("^CON_", isotab$protein ) , ]
isotab <- isotab[ !grepl("^REV_", isotab$protein ) , ]
isotab <- isotab[ !grepl("^NA-", isotab$protein ) , ]
isotab$nondis <- samtot-isotab$samples

isotab$protein <- factor(isotab$protein, levels = isotab$protein[order(isotab$samples)])
isotab <- melt(isotab)


ggplot(isotab, aes(x = protein, y = value, fill=variable)) +
  geom_bar(stat='identity', position = position_fill(reverse = TRUE) ) +
  coord_flip() +
  ylab("fraction of samples") +
  scale_fill_manual("", values = c("darkgreen", "grey"), labels = c("detected", "not detected")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
  theme(legend.position="bottom")

width=18.5
ggsave("./figures/proteins_recovered_ISOFORMS_proteinGroups.pdf", width = width, height = 3.5 * width, limitsize = FALSE)

