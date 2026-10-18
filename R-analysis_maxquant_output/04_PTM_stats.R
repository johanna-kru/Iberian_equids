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


######
### Read in data
######

  setwd("Ranalysis_RUN9")  

  ### Ages
  ages <- read.csv(file = "samples_ages.csv", sep ="\t" , stringsAsFactors = FALSE)
  
######
### Deamidation separate because it needs combination of N and Q
######

### Files conversion rates

  ### Function to create bbr_df
  bbr_create <- function(ptm){ 
    
    fnb <- read.table(paste("./evidence_all/" , ptm ,"/filenames_bulk_bootstrap_replicates.txt" ,sep = ""), stringsAsFactors = FALSE)$V1
    
    
  bbr <- list()
  bbr_df <- data.frame()
  ### bulk bootstrap
    for(ii in fnb){
      
      bbr[[ii]] <- read.csv(paste( "./evidence_all/" ,ptm,"/", ii, sep = "") , stringsAsFactors = FALSE)
    
      # get sample name
      st <- bbr[[ii]]$Raw_file[1]
      
      # Remove unwanted characters
      st <- gsub("[(',)]", "", st)
      
      # Split the string by underscore
      split_st <- strsplit(st, "_")
      
      # Get the last element
      last_element <- tail(split_st[[1]], n = 1)
      
      bbr[[ii]]$samples <- last_element
      
      bbr_df <- rbind(bbr_df, bbr[[ii]][, c( "Amino_acid", "Conversion_ratio", "samples"  ) ])
      
          }

  bbr$bulk_bootstrap_replicates_P1327B_FT.csv$samples <- "P1327B_FT"
  
  bbr_df <- merge( bbr_df, ages, by="samples", all.y = TRUE)
  bbr_df$samples <- factor(bbr_df$samples, levels = unique( bbr_df$samples[order(bbr_df$ageBP_Ma)]) )
  
  return(bbr_df)
  }
  
 bbr_dfQ <- bbr_create("deamidationQ")
 bbr_dfN <- bbr_create("deamidationN")
 
 ### include Arginine to Ornintine
 bbr_dfR <- bbr_create("arg2orn")

 ### include Arginine to Ornintine  
 bbr_df <- rbind(bbr_dfQ, bbr_dfN,bbr_dfR)
 
 
  ### Function to read per protein conversion
 
 ppc_create <- function(ptm){
   
  fnp <- read.table(paste("./evidence_all/" , ptm ,"/filenames_per_protein_conversion_ratio.txt" ,sep = ""), stringsAsFactors = FALSE)$V1

  ppc <- list()
  ppcagg <- list()
  
  for(ii in fnp){
    ppc[[ii]] <- read.csv(paste( "./evidence_all/" ,ptm,"/", ii, sep = "")  , stringsAsFactors = FALSE)
    
    # get sample name
    st <- ppc[[ii]]$Raw_file[1]
    
    # Remove unwanted characters
    st <- gsub("[(',)]", "", st)
    
    # Split the string by underscore
    split_st <- strsplit(st, "_")
    
    # Get the last element
    last_element <- tail(split_st[[1]], n = 1)
    
    ppc[[ii]]$samples <- last_element
    
    ppc[[ii]]$Leading_razor_protein <- sapply(ppc[[ii]]$Leading_razor_protein, function(x) sub("^SAV_", "", x))
    
    pp <- strsplit( ppc[[ii]]$Leading_razor_protein , "_" )
    
    pp <- sapply(pp, function(x) x[1])
    
    ppc[[ii]]$protein <- pp
    
  }
  
  ppc$per_protein_conversion_ratio_P1327B_FT.csv$samples <- "P1327B_FT"

  
  ### aggregate proteins
  for(ii in fnp){
     ppcagg[[ii]] <- ppc[[ii]] %>%
      select(-Raw_file, -Leading_razor_protein) %>%
      group_by(protein) %>%
      summarise(
        samples = first(samples),
        PTM = first(PTM),
        Amino_acid_count = sum(Amino_acid_count),
        Conversion_count = sum(Conversion_count),
        Abundance_count = sum(Abundance_count),
        Conversion_ratio = mean(Conversion_ratio),
        Boostrap_mean_converion_ratio = mean(Boostrap_mean_converion_ratio),
        Q95_low = mean(Q95_low),
        Q95_up = mean(Q95_up)
      )

       }
  
  ppcagg_df <- do.call(rbind, ppcagg)
  
  return(ppcagg_df)
  
 }
  
 ppcagg_dfQ <- ppc_create("deamidationQ")
 ppcagg_dfN <- ppc_create("deamidationN")
 
 ### include Arginine to Ornintine
 ppcagg_dfR <- ppc_create("arg2orn")
 
 ppcagg_dfQ$Amino_acid <- "Q"
 ppcagg_dfN$Amino_acid <- "N"
 ### include Arginine to Ornintine
 ppcagg_dfR$Amino_acid <- "R"
 ### include Arginine to Ornintine
 ppcagg_df <- rbind(ppcagg_dfQ, ppcagg_dfN, ppcagg_dfR)
 
 obs <- aggregate(Amino_acid_count ~ samples + Amino_acid, ppcagg_df, sum, na.rm = TRUE)


 ### plot bulk bootstrap
  
dir.create("figures")
  
  # subset to samples produced in this study
  firstsam <-  bbr_df[ grepl("P.*A", bbr_df$samples , perl = TRUE )  , "samples"]
  
  # exchange CPH sample (P1327A) for BCN version of same individual (P1327B)
  firstsam[firstsam == "P1327A"] <- "P1327B"

  # E.hydruntinus specimens (teeth) from same individual, compare only one sample from one tooth per individual
  firstsam<- firstsam[ !( firstsam %in% c("P2856A", "P2857A") ) ]
  
  bbr_df <- bbr_df[bbr_df$samples %in% firstsam, ]
  lvs <- unique( bbr_df$samples[order(bbr_df$ageBP_Ma)])
  lvs <- lvs[!lvs %in% "P3006A"]
  bbr_df$samples <- factor(bbr_df$samples, levels = lvs )
  obs <- obs[obs$samples %in% firstsam,]
  
    # Merge the dataframes
  merged_df <- merge(bbr_df, obs, by = c("samples", "Amino_acid"))

  # Ensure the levels of the 'samples' factor in 'obs' match those in 'bbr_df'
  obs$samples <- factor(obs$samples, levels = levels(merged_df$samples))
  # Plot
  ggplot(data = merged_df, aes(x = samples, y = Conversion_ratio, fill = Amino_acid)) +
    geom_violin(aes(color = Amino_acid), position = position_dodge(width = 0.9)) +
    geom_text(data = obs, aes(label = Amino_acid_count, y = 1.1), size = 3 * 1.3, angle = 0, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = c( mpal[1], mpal[2], mpal[3] )) +  # Choose your own colors  #include Arginine to Ornintine
    scale_color_manual(values = c(mpal[1], mpal[2], mpal[3] )) +  # Outline colors  # include Arginine to Ornintine
    ylab("Conversion Ratio") +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "lightgrey"),
          panel.border = element_rect(colour = "lightgrey", fill=NA, linewidth=1.4))+
    scale_y_continuous(
      limits = c(-0.05, 1.2),
      expand = expansion(mult = c(0, 0.05))
    )+
    coord_flip() +
    scale_x_discrete(limits = rev(levels(merged_df$samples)))
    

  
  
  width = 9
   ### include Arginine to Ornintine
  ggsave("./figures/bulk_bootstrap_deamidation_arg2orn.pdf", width = 0.5 * width, height = width)
  
      

  ### plot by protein
  
  # subset to samples produced in this study
  ppcagg_df <- merge( ppcagg_df, ages, by="samples", all.y = TRUE)
  ppcagg_df <- ppcagg_df[ppcagg_df$samples %in% firstsam, ]
  ppcagg_df[ppcagg_df$samples == "P3006A", "protein"] <- "AMELX" # set protein to some value to avoid creation of NA as another categors in proteins
  ppcagg_df$Amino_acid <- factor(ppcagg_df$Amino_acid, levels = c("Q", "N")) # needed to dodge the two bars in plot next top each other
  ppcagg_df$samples <- factor(ppcagg_df$samples, levels = unique( ppcagg_df$samples[order(ppcagg_df$ageBP_Ma)]) )
    
  ggplot(data = ppcagg_df, aes(x = protein, y = Conversion_ratio , fill = Amino_acid)  )+
    geom_bar(stat='identity', position = position_dodge2(preserve = "single"), width=0.7) + #preserve = "single" so that single barplot not size of both
    scale_fill_manual(values = alpha( c( mpal[2], mpal[3] ), 0.8)) +  # Choose your own colors
    geom_errorbar(aes(ymin=Q95_low, ymax=Q95_up), width=0.5, linewidth=0.5, color="grey28", position = position_dodge2(preserve = "single")) + # also dodge errorbars
    geom_text(aes(label = Amino_acid_count, y = 0.1), size = 2.7, angle = 90, position = position_dodge(width = 0.9), colour="black") +
    facet_wrap(samples ~ ., ncol = 3, strip.position = "right") +
    ylab("Conversion Ratio") +
    theme(axis.text.x = element_text( angle = 45, hjust = 1) ) +
    theme(axis.text.x = element_text(size = 10)) +
    theme( panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), 
           panel.background = element_blank(),axis.line = element_line(colour = "lightgrey"),
          panel.border = element_rect(colour = "lightgrey", fill=NA, linewidth=1.4)) +
    scale_y_continuous(n.breaks = 6 ) + 
    theme(legend.position = "bottom")

  
  
  width = 12
  ggsave("./figures/per_protein_conversion_deamidation.pdf", width = width, height = 0.9 * width)
  
  
  
  ######
  ### Function to plot any conversion ratio
  ######
  
  
  ### Files conversion rates
  
  conrat <- function(ptm){

  
  fnb <- read.table(paste("./evidence_all/" , ptm ,"/filenames_bulk_bootstrap_replicates.txt" ,sep = ""), stringsAsFactors = FALSE)$V1
  
  fnp <- read.table(paste("./evidence_all/" , ptm ,"/filenames_per_protein_conversion_ratio.txt" ,sep = ""), stringsAsFactors = FALSE)$V1
  
  bbr <- list()
  bbr_df <- data.frame()
  ppc <- list()
  ppcagg <- list()
  proteins <- list()
  
  ### bulk bootstrap
  for(ii in fnb){
    
    bbr[[ii]] <- read.csv(paste( "./evidence_all/" ,ptm,"/", ii, sep = ""), stringsAsFactors = FALSE)
    
    # get sample name
    st <- bbr[[ii]]$Raw_file[1]
    
    # Remove unwanted characters
    st <- gsub("[(',)]", "", st)
    
    # Split the string by underscore
    split_st <- strsplit(st, "_")
    
    # Get the last element
    last_element <- tail(split_st[[1]], n = 1)
    
    bbr[[ii]]$samples <- last_element
    
    bbr_df <- rbind(bbr_df, bbr[[ii]][, c( "Amino_acid", "Conversion_ratio", "samples"  ) ])
    
  }
  
  bbr$bulk_bootstrap_replicates_P1327B_FT.csv$samples <- "P1327B_FT"
  
  bbr_df <- merge( bbr_df, ages, by="samples", all.y = TRUE)
  bbr_df$samples <- factor(bbr_df$samples, levels = unique( bbr_df$samples[order(bbr_df$ageBP_Ma)]) )
  
  
  ### per protein conversion
  for(ii in fnp){
    ppc[[ii]] <- read.csv(paste( "./evidence_all/" ,ptm,"/", ii, sep = "") , stringsAsFactors = FALSE)
    
    # get sample name
    st <- ppc[[ii]]$Raw_file[1]
    
    # Remove unwanted characters
    st <- gsub("[(',)]", "", st)
    
    # Split the string by underscore
    split_st <- strsplit(st, "_")
    
    # Get the last element
    last_element <- tail(split_st[[1]], n = 1)
    
    ppc[[ii]]$samples <- last_element
    
    ppc[[ii]]$Leading_razor_protein <- sapply(ppc[[ii]]$Leading_razor_protein, function(x) sub("^SAV_", "", x))
    
    pp <- strsplit( ppc[[ii]]$Leading_razor_protein , "_" )
    
    pp <- sapply(pp, function(x) x[1])
    
    ppc[[ii]]$protein <- pp
    
    proteins[[ii]] <- unique(sort(pp))
    
  }
  
  ppc$per_protein_conversion_ratio_P1327B_FT.csv$samples <- "P1327B_FT"
  
  
  ### aggregate proteins
  for(ii in fnp){
    ppcagg[[ii]] <- ppc[[ii]] %>%
      select(-Raw_file, -Leading_razor_protein) %>%
      group_by(protein) %>%
      summarise(
        samples = first(samples),
        PTM = first(PTM),
        Amino_acid_count = sum(Amino_acid_count),
        Conversion_count = sum(Conversion_count),
        Abundance_count = sum(Abundance_count),
        Conversion_ratio = mean(Conversion_ratio),
        Boostrap_mean_converion_ratio = mean(Boostrap_mean_converion_ratio),
        Q95_low = mean(Q95_low),
        Q95_up = mean(Q95_up)
      )
    
  }
  
  ppcagg_df <- do.call(rbind, ppcagg)
  
  obs <- aggregate(Amino_acid_count ~ samples, ppcagg_df, sum, na.rm = TRUE)
  
  
  ### plot bulk bootstrap
  
  dir.create("figures")
  
  #  bbr_df <- do.call(rbind, bbr)
  
  # subset to samples produced in this study
  firstsam <-  bbr_df[ grepl("P.*A", bbr_df$samples , perl = TRUE )  , "samples"]
  
  # exchange CPH sample (P1327A) for BCN version of same individual (P1327B)
  firstsam[firstsam == "P1327A"] <- "P1327B"
  # E.hydruntinus specimens (teeth) from same individual, compare only one sample from one tooth per individual
  firstsam<- firstsam[ !( firstsam %in% c("P2856A", "P2857A") ) ]
  
  bbr_df <- bbr_df[bbr_df$samples %in% firstsam, ]
  bbr_df$samples <- factor(bbr_df$samples, levels = unique( bbr_df$samples[order(bbr_df$ageBP_Ma)]) )
  obs <- obs[obs$samples %in% firstsam,]
  
  # Merge the dataframes
  merged_df <- merge(bbr_df, obs, by = c("samples"))
  lvs <- unique( bbr_df$samples[order(bbr_df$ageBP_Ma)])
  lvs <- lvs[!lvs %in% "P3006A"]
  merged_df$samples <- factor(merged_df$samples, levels =  lvs)
  # Ensure the levels of the 'samples' factor in 'obs' match those in 'bbr_df'
  obs$samples <- factor(obs$samples, levels = levels(merged_df$samples))
  
  ggplot(data = merged_df, aes(x = samples, y = Conversion_ratio) )+
    geom_violin( fill = mpal[2], color = mpal[2]) +
    #    geom_boxplot(width = 0.01, color = "darkgrey", outliers = FALSE) + 
    geom_text(data = obs, aes(label = Amino_acid_count, y = -0.02), size = 3 * 1.3, angle = 0, position = position_dodge(width = 0.9)) +
    ylab("Conversion Ratio") +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_blank() )+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "lightgrey"),
          panel.border = element_rect(colour = "lightgrey", fill=NA, linewidth=1.4)) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(merged_df$samples)))
  
  width = 9
  ggsave( paste( "./figures/", "bulk_bootstrap_" , ptm , ".pdf", sep=""), width = 0.5 * width, height = width)
  
  
  
  ### plot by protein
  
  # subset to samples produced in this study
  ppcagg_df <- merge( ppcagg_df, ages, by="samples", all.y = TRUE)
  ppcagg_df <- ppcagg_df[ppcagg_df$samples %in% firstsam, ]
  ppcagg_df$samples <- factor(ppcagg_df$samples, levels = unique( ppcagg_df$samples[order(ppcagg_df$ageBP_Ma)]) )

  
  ggplot(data = ppcagg_df, aes(x = protein, y = Conversion_ratio) )+
    geom_bar(stat='identity', position= position_dodge(0.7), width=0.7, fill = mpal[2], alpha = 0.8) +
    geom_errorbar(aes(ymin=Q95_low, ymax=Q95_up), width=0.5, linewidth=0.5, color="grey28") +
    geom_text(aes(label = Amino_acid_count, y = 0.1), size = 2.7, angle = 0, position = position_dodge(width = 0.9), colour="black") +
    facet_wrap(samples ~ ., ncol = 2, strip.position = "right") +
    ylab("Conversion Ratio") +
    theme(axis.text.x = element_text( angle = 45, hjust = 1) ) +
    theme(axis.text.x = element_text(size = 10)) +
    theme( panel.grid.major = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), 
           panel.background = element_blank(),axis.line = element_line(colour = "lightgrey"),
           panel.border = element_rect(colour = "lightgrey", fill=NA, linewidth=1.4)) +
    scale_y_continuous(n.breaks = 6 ) 
  
  
  
  width = 9
  ggsave( paste( "./figures/", "per_protein_conversion_" , ptm , ".pdf", sep=""), width = width, height = 1.2 * width)
  
  
  }
  
  conrat("arg2orn")
  conrat("Gln2pyro-Glu")
  conrat("Glu2pyro-Glu")
  
  