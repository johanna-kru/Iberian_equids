rm(list=ls())

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggridges)
library(RColorBrewer) # colors 
library("viridis") # colors 



### define colors for most plots

### master colors
mpal <- c("#2853a8ff","#7E97CBFF", 
          "#c32f31ff", "#FB8072" , 
          "#2f8450ff", "#62c1b8ff" , 
          "#7b43a9ff","#b98ed7ff", 
          "chocolate3","#FDB462" )

#338CC8

### colors for plot in supplement
pal <- c("#FDBF6F", "#E31A1C", "#FB9A99", "darkgreen", "#A6CEE3" , "#1F78B4",  "mediumpurple", "#FFED6F", "#CCEBC5", "#BC80BD", "#FCCDE5",
         "#B3DE69" , "#FDB462", "#80B1D3", "#FB8072", "#BEBADA", "#FFFFB3", "#8DD3C7", "#D9D9D9" )

######
### Read in data
######

setwd("Ranalysis_RUN9")

### Ages
ages <- read.csv(file = "samples_ages.csv", sep ="\t" , stringsAsFactors = FALSE)


### Link between IDs and gene names
id2gn <- read.table(file = "id2genename_SAP_test_db_BCN_horses.txt", stringsAsFactors = FALSE)
colnames(id2gn) <- c( "protein", "id")
id2gn <- id2gn[!duplicated(id2gn), ]


### Amino acid counts SPECIMENS WITH SEVERAL SAMPLES CONCATENATED

aacounts <- read.table("./consensus/concatenated_by_specimen/aa_counts.tsv", sep = '\t', stringsAsFactors = FALSE)
colnames(aacounts) <- ( c("samples", "protein", "aa") )
aacounts <- aacounts[aacounts$protein != "SPI2",]


### Coverage
  fl <-  read.table(file = "./coverage/concatenated_by_specimen/filenames.txt", stringsAsFactors = FALSE)$V1
  
  covl <- list()
  
  for(ii in fl){
    
    print(ii)
    fstem <-  sub(".tsv", "", sub("coverage_", "", ii) ) 
    protein <- sample <- strsplit(fstem, "_")[[1]][1]
    sample <- strsplit(fstem, "_")[[1]][2]
 
    # Read the file
    file_content <- readLines( paste( "./coverage/concatenated_by_specimen/", ii, sep = "") )
    
    # Check if the file is empty
      if(length(file_content) == 0) {
        addf <-data.frame(position=NA, coverage=NA)
      } else {
  
      addf <- try( read.table( paste( "./coverage/concatenated_by_specimen/", ii, sep = "") , header = TRUE, sep = '\t', stringsAsFactors = FALSE) )
      }  


    
    colnames(addf)[colnames(addf) == "coverage"] <- paste("coverage_", sample, sep = "")
    
    covl[[ii]] <- addf
      
  }
  
  
######  
### Merge coverage by protein
######
  
  # Get the names of the list elements
  names <- names(covl)
  
  # Extract the categories from the names
  categories <- gsub(".*_(.*)_.*", "\\1", names)

  # Get the unique categories
  proteins <- unique(categories)
  
  # Initialize an empty list to store the merged data frames
  merged_cov <- list()
  cov_all <- data.frame()
  
    # For each unique protein
    for (protein in proteins) {
      # Get the names of the list elements that belong to the current protein
      protein_names <- names[categories == protein]
      
      # Get the data frames that belong to the current protein
      protein_dfs <- covl[protein_names]
      max_rows <- max(sapply(protein_dfs, nrow))
      
      covl_padded <- lapply(protein_dfs, function(df) {
        num_rows <- nrow(df)
        if (num_rows < max_rows) {
          # Create a data frame of NAs with the correct number of rows and columns
          df_na <- as.data.frame(matrix(NA, nrow = max_rows - num_rows, ncol = ncol(df)))
          colnames(df_na) <- colnames(df)
          
          # Add the rows of NAs to the data frame
          df <- rbind(df, df_na)
          df$position <- c(1:max_rows)
        }
        return(df)
      })
      
      # Merge the data frames by the "position" column
      merged_df <- Reduce(function(x, y) merge(x, y, by = "position"), covl_padded)
      
      merged_df$sum_cov <- rowSums(merged_df[,!(colnames(merged_df) %in% c("position") ) ], na.rm = TRUE)
      
      merged_df$protein <- protein
      
      # Store the merged data frame in the list
      merged_cov[[protein]] <- merged_df
      
      cov_all <- rbind(cov_all, merged_df[, c("position", "sum_cov", "protein")] )
    }
    
    cov_all <- cov_all[ is.na(cov_all$position) == FALSE,]
    cov_all[cov_all == 0] <- NA
  
    
    # independent x and y scales
    ggplot(cov_all, aes(x = position, y = sum_cov)) +
      geom_bar(stat='identity', position= position_dodge(0), width=1, fill = mpal[1]) +
      ylab("number PSMs") +
      facet_wrap(protein ~ ., scales="free", ncol = 3, strip.position = "right") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "lightgrey"),
            panel.border = element_rect(colour = "lightgrey", fill=NA, linewidth=1.4)) +
      theme(axis.title = element_text(size = 13)) +
      theme(axis.text.x = element_text(size = 10)) +
      theme(axis.text.y = element_text(size = 10)) +
      
      theme(legend.position="none")
    
    
    width = 7
    ggsave("./figures/coverage_by_protein_CONCATENATED_BY_SPECIMEN.pdf", width = width, height = 0.8 * width)
    

    
    # log10 y-scale
  ggplot(cov_all, aes(x = position, y = sum_cov)) +
    geom_bar(stat='identity', position= position_dodge(0), width=1, fill = mpal[1]) +
    scale_y_log10() +
    ylab("number PSMs") +
    facet_wrap(protein ~ ., scales="free", ncol = 3, strip.position = "right") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "lightgrey"),
          panel.border = element_rect(colour = "lightgrey", fill=NA, linewidth=1.4)) +
    theme(axis.title = element_text(size = 13)) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 10)) +

    theme(legend.position="none")
  
  
  width = 7
  ggsave("./figures/coverage_by_protein_log10_CONCATENATED_BY_SPECIMEN.pdf", width = width, height = 0.8 * width)
  
  

  
  # consistent y-scale
  ggplot(cov_all, aes(x = position, y = sum_cov)) +
    geom_bar(stat='identity', position= position_dodge(0), width=1, fill = mpal[1]) +
    ylab("number PSMs") +
    facet_wrap(protein ~ ., scales="free_x", ncol = 3, strip.position = "right") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "lightgrey"),
          panel.border = element_rect(colour = "lightgrey", fill=NA, linewidth=1.4)) +
    theme(axis.title = element_text(size = 13)) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 10)) +
    
    theme(legend.position="none")
  
  
  width = 7
  ggsave("./figures/coverage_by_protein_consistent_y_scale_CONCATENATED_BY_SPECIMEN.pdf", width = width, height = 0.8 * width)

  
  
  
######
### Plot total amount of amino acids by sample ALL PROTEINS
######
  
  # merge with ages (only account for first 6 characters)
  
  # Create new columns with the first 6 characters of 'samples'
  aacounts$merge_id <- substr(aacounts$samples, 1, 6)
  aacounts$merge_id[ aacounts$merge_id == "P2855-"] <- "P2855A"
  ages$merge_id     <- substr(ages$samples, 1, 6)
  
  # Merge using the new column and rename
  aacounts <- merge(aacounts, ages, by = "merge_id", all = FALSE)
  aacounts <- aacounts[, c("samples.x", "protein", "aa", "ageBP_Ma") ]
  colnames(aacounts )[colnames(aacounts ) == "samples.x"] <- "samples"

  aacounts$samples <- factor(aacounts$samples, levels = unique( aacounts$samples[order(aacounts$ageBP_Ma)]) )
  
  
  # subset to samples produced in this study (in this case all specimens)
  firstsam <-  aacounts[ grepl("P.*", aacounts$samples , perl = TRUE )  , "samples"]
  
  # exchange CPH sample (P1327A) for BCN version of same individual (P1327B)
  if (!"P1327B" %in% firstsam) {
    firstsam <- c(as.character(firstsam), "P1327B")
  }
  
  # E.hydruntinus specimens (teeth) from same individual, compare only one sample from one tooth per individual
  firstsam<- firstsam[ !( firstsam %in% c("P2856A", "P2857A") ) ] 
  
  
  totalaaFirst <- aacounts[ aacounts$samples %in% firstsam  ,]
  totalaaFirst$samples <- factor(totalaaFirst$samples, levels = unique( totalaaFirst$samples[order(totalaaFirst$ageBP_Ma)]) )
  totalaaFirst_agg <- aggregate(aa ~ samples, data=totalaaFirst ,FUN=sum) #aggregation manual, otherwise counts divided by ten in ggplot
  
  ggplot(totalaaFirst_agg, aes(x = samples, y = aa)) +
    geom_bar(stat='identity', position= position_dodge(0.7), width=0.7, fill = mpal[2]) +
    ylab("amino acids recovered") +
    scale_x_discrete(limits = rev(levels(totalaaFirst$samples))) + # needed to reverse order of factors
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
    theme(axis.title = element_text(size = 28)) +
    theme(axis.text.x = element_text(size = 29)) +
    theme(axis.text.y = element_text(size = 29)) +
    theme(legend.position="none")
  
  
  width = 12
  ggsave("./figures/total_aa_ASRs_CONCATENATED_BY_SPECIMEN.pdf", width = width, height = 1.6 * width)
  
  
  ######
  ### Relative frequencies aa by protein ALL PROTEINS
  ######
  
  ggplot(totalaaFirst, aes(x = samples, y = aa, fill = protein)) +
    geom_bar(stat='identity', position="fill", width=0.7, ) +
    scale_x_discrete(limits = rev(levels(totalaaFirst$samples))) + # needed to reverse order of factors
    coord_flip() +
    scale_fill_manual(values = rep(pal, length.out=19)) +
    scale_y_reverse() +
    ylab("fraction amino acids by protein") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
    theme(axis.title = element_text(size = 15)) +
    theme(axis.text.x = element_text(size = 13)) +
    theme(axis.text.y = element_text(size = 13)) 
  
  
  width = 12
  ggsave("./figures/rel_ACCEPTED_aa_CONCATENATED_BY_SPECIMEN.pdf", width = width, height = 0.6 * width)
  
  
  
  ######
  ### Plot total amount of amino acids by sample NO COLLAGENS ANCIENT
  ######
  
  collagens <- c("COL1A1", "COL1A2", "COL2A1")
  samex <- c("P2859A", "P2860A", "P2861A")
  totalaaFirstnoCol <- totalaaFirst[!( (totalaaFirst$samples %in% samex) & (totalaaFirst$protein %in% collagens) ),]
  totalaaFirstnoCol_agg <- aggregate(aa ~ samples, data=totalaaFirstnoCol ,FUN=sum) #aggregation manual, otherwise counts divided by ten in ggplot
  
  ggplot(totalaaFirstnoCol_agg, aes(x = samples, y = aa)) +
    geom_bar(stat='identity', position= position_dodge(0.7), width=0.7, fill = mpal[2]) +
#   scale_fill_manual(values = bluehue) +
    ylab("amino acids recovered") +
    scale_x_discrete(limits = rev(levels(totalaaFirstnoCol$samples))) + # needed to reverse order of factors
    coord_flip() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
    theme(axis.title = element_text(size = 28)) +
    theme(axis.title.y = element_blank()) +
    theme(axis.text.x = element_text(size = 29)) +
    theme(axis.text.y = element_text(size = 29)) +
    theme(legend.position="none")
  
  
  width = 14
  ggsave("./figures/total_aa_ASRs_NO_COLLAGENS_CONCATENATED_BY_SPECIMEN.pdf", width = width, height = 19.2)
  
  
  ######
  ### Relative frequencies aa by protein ALL PROTEINS NO COLLAGENS ANCIENT
  ######
  
  ggplot(totalaaFirstnoCol, aes(x = samples, y = aa, fill = protein)) +
    geom_bar(stat='identity', position="fill", width=0.7, ) +
    scale_x_discrete(limits = rev(levels(totalaaFirstnoCol$samples))) + # needed to reverse order of factors
    coord_flip() +
    scale_fill_manual(values = rep(pal, length.out=19)) +
    scale_y_reverse() +
    ylab("fraction amino acids by protein") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
    theme(axis.title = element_text(size = 15)) +
    theme(axis.text.x = element_text(size = 13)) +
    theme(axis.text.y = element_text(size = 13)) 
  
  
  width = 12
  ggsave("./figures/rel_ACCEPTED_aa_NO_COLLAGENS_CONCATENATED_BY_SPECIMEN.pdf", width = width, height = 0.6 * width)
  
  ######
  ### Relative frequencies aa by selected proteins ALL PROTEINS NO COLLAGENS ANCIENT
  ######
  
### combine less frequent proteins to "other"
otherprot <- c("AHSG", "ALB", "AMTN", "COL2A1", "KLK4", "MGP", "ODAM", "SERPINA1", "SERPINC1", "TTR")
  
totalaaFirstnoCol_compact <- totalaaFirstnoCol

totalaaFirstnoCol_compact[totalaaFirstnoCol_compact$protein %in% otherprot, "protein"] <- "other" 

### merge all entries with "other"
df_other <- totalaaFirstnoCol_compact %>%
  filter(protein == "other") %>%
  group_by(samples) %>%
  summarise(protein = "other", aa = sum(aa), ageBP_Ma = first(ageBP_Ma))

df_main <- totalaaFirstnoCol_compact %>%
  filter(protein != "other")

totalaaFirstnoCol_compact <- bind_rows(df_main, df_other) %>%
  arrange(samples, desc(protein))



ggplot(totalaaFirstnoCol_compact, aes(x = samples, y = aa, fill = protein)) +
  geom_bar(stat='identity', position="fill", width=0.7, ) +
  scale_x_discrete(limits = rev(levels(totalaaFirstnoCol$samples))) + # needed to reverse order of factors
  coord_flip() +
  scale_fill_manual(values = c(mpal[9], mpal[1:7],"lightgrey" ) ) +
  scale_y_reverse() +
  ylab("fraction amino acids by protein") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "lightgrey")) +
  theme(axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 13)) +
  theme(axis.text.y = element_text(size = 13)) 


width = 12
ggsave("./figures/rel_ACCEPTED_aa_NO_COLLAGENS_compact_CONCATENATED_BY_SPECIMEN.pdf", width = width, height = 0.6 * width)

  