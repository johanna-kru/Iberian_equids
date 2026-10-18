

######
### Read in data
######

### list of filenames
fm <- read.table(file = "./msms_all/filenames.txt", stringsAsFactors = FALSE)$V1

### define samples
samples <-  sub(".txt", "", sub("msms_", "", fm) ) 

### link between IDs and gene names
id2gn <- read.table(file = "id2genename_SAP_test_db_BCN_horses.txt", stringsAsFactors = FALSE)
#id2gn <- id2gn[, c(1,7)] #for older versions 
colnames(id2gn) <- c( "protein", "id")
id2gn <- id2gn[!duplicated(id2gn), ]

### read in msms files
msms <- list()
msms_sel <- list()

for(ii in fm){
  
  print(ii)
  sample <-  sub(".txt", "", sub("msms_", "", ii) ) 
  
  ### define relevant columns, spaces read in a s "." in R
  relsel <- c( "Scan.number", "Scan.index", "id", "Peptide.ID", "Sequence", "Length", "Proteins" , "mainID", "PEP", "Score", "Delta score", "Matches")
  
  
  # use try so that loop continues if one file empty
  msms[[ii]] <- try (read.csv( paste( "./msms_all/", ii, sep = "")  , sep = '\t', stringsAsFactors = FALSE) )
  
  # sapply function to predefined column of df
  msms[[ii]]$mainID <- sapply( msms[[ii]]$Proteins, function(x)unlist(strsplit( as.character(x), ";" )) [1] )
  
  
  # subset to important columns
  msms_sel[[ii]] <- msms[[ii]][, names(msms[[ii]]) %in% as.character( relsel ) ]
  #  msms_sel[[ii]]$No.Scans <- sapply( msms_sel[[ii]]$MS.MS.IDs , function(x)length( strsplit(as.character(x), split = ";")[[1]] ) )
  
  msms_sel[[ii]][msms_sel[[ii]] == 0] <- NA
  
  ### Link IDs to protein names (without species, if species wanted add columns 3,4 in id2gn). 
  ###REV and CON entries won't have a "protein" entry, so that "protein.long" will start with NA. This is needed to eliminate them later!
  msms_sel[[ii]] <- merge(x = msms_sel[[ii]], y = id2gn[, c(1,2)], all.x = TRUE, by.x = "mainID", by.y = "id")
  msms_sel[[ii]]$protein.long <- paste( msms_sel[[ii]]$protein, msms_sel[[ii]]$mainID, sep = "-" ) 
  
  ### Replace missing gene names with Protein.IDs (for contaminants)
  msms_sel[[ii]]$protein[ is.na( msms_sel[[ii]]$protein ) ] <- msms_sel[[ii]]$mainID[ is.na( msms_sel[[ii]]$protein ) ]

  ### Remove rows where Matches is empty
  msms_sel[[ii]] <- msms_sel[[ii]][msms_sel[[ii]]$Matches != "", ]
  
}

#View(msms_sel$msms_P1327B.txt)



######
### Ion series completeness score
######

min_peplength <- 7
# peptide of 7aa has 6 bonds if full coverage isc= 6/6 0 1.0
# If one missing in shortest possible peptide isc = 5/6 = 0.833 -> define as minimum here

for( ii in fm ){
  
  
  icScore <- c()
  bionScore <- c()
  yionScore <- c()
  ricScore <- c()
  bonds <- c()
  remark <- c()
  
  for( jj in 1:dim(msms_sel[[ii]])[1] ){
    
    ll <- msms_sel[[ii]]$Length[jj]
    
    sr <- msms_sel[[ii]]$Matches[jj]
    sr <- strsplit(sr, ";")[[1]]
    sr <- sub( "-.*" , "" , sr, perl=TRUE)
    sr <- sub("\\*", "", sr)
    sr <- sub("\\(.*", "", sr)
    rma <- grepl("^a", sr)
    sr <- sr[!rma]
    sr <- unique( sort(sr) )
    
    dfsr <- as.data.frame(sr)
    dfsr$type <- gsub("[[:digit:]]","",dfsr$sr)
    dfsr$index <- gsub("[^0-9.-]", "",dfsr$sr)
    dfsr$pepbond <- as.numeric(dfsr$index)
    
    ### 1. raw ic score
    
    ### convert everything as if b ions
    
    rr <- dfsr$type != "b"
    converted_index <- rep(ll, length(rr[rr == TRUE]) ) - as.numeric( dfsr[rr, "index"] )
    
    dfsr[ rr , "pepbond"  ] <- converted_index
    dfsr$converted_ions <- paste( "i", dfsr$pepbond, sep = "")
    
    conv_ions <-  unique( sort(dfsr$converted_ions) )
    
    ### subtract -1 because max no of peptide bonds = residues -1
    ic <- format( round( length(conv_ions)/(ll -1) , digits = 4), nsmall = 4) 
    
    icScore <- c(icScore, ic)
    
    
    
    ### 2. b-ion score
    
    b_ions <- dfsr[dfsr$type == "b", "sr"]
    b_ions <- unique( sort(b_ions) )
    bc <- format( round( length(b_ions)/(ll -1) , digits = 4), nsmall = 4) 
    bionScore <- c(bionScore, bc)
    
    
    ### 3. y-ion score
    y_ions <- dfsr[dfsr$type == "y", "sr"]
    y_ions <- unique( sort(y_ions) )
    yc <- format( round( length(y_ions)/(ll -1) , digits = 4), nsmall = 4) 
    yionScore <- c(yionScore, yc)
    
    
    ### 4. relaxed score that forgives if b(1), b(ll-1) and y(1), y(ll-1) missing but penalty if missing in between
    
    # Check if b(1) or y(ll-1) exist
    if( any( dfsr$pepbond == (ll-1) ) == TRUE ){
      
      pl <- ll 
      remrk <- " b(1) or y(n-1) exist. "
    } else {
      
      pl <- ll-1 
      remrk <- " NO b(1) or y(n-1). "
      
    }
    
    
    
    # y(1) or b(ll-1) exist
    if( any( dfsr$pepbond == 1 ) == TRUE ){
      
      ml <- pl 
      
      remrk <- paste(remrk, " y(1) or b(n-1) exist. ", sep = "")
      
    }else{
      
      ml <- pl -1 
      remrk <- paste(remrk, " NO y(1) or b(n-1). ", sep = "")
      
    }
    
    remark <- c(remark, remrk)
    
    
    ric <- format( round( length(conv_ions)/(ml-1) , digits = 4), nsmall = 4) 
    ricScore <- c(ricScore, ric)
    
    ### if shortest (relaxed) peptide possible (in most cases 5 aa), one missing peptide bond results in ic = 3/4 = 0.75
    
    
    ### 5. Create string that indicates how well covered each peptide bond is
    # Create an empty character vector of length maxpepbond
    
    maxpepbond <- ll-1
    output <- rep("X", maxpepbond)
    
    # For each unique value in pepbond
    for (i in unique(dfsr$pepbond)) {
      # Subset the dataframe for rows where pepbond equals i
      subset_dfsr <- dfsr[dfsr$pepbond == i,]
      
      # Check if there is a "b" or a "y" in the "type" column
      if ("b" %in% subset_dfsr$type & "y" %in% subset_dfsr$type) {
        output[i] <- "P"
      } else if ("b" %in% subset_dfsr$type) {
        output[i] <- "b"
      } else if ("y" %in% subset_dfsr$type) {
        output[i] <- "y"
      }
    }
    
    # Convert the output vector to a single string
    output_string <- paste(output, collapse = "")
    
    # Print the output string
    bonds <- c(bonds, output_string)
    
    
  }
  
  msms_sel[[ii]]$icScore <- icScore
  msms_sel[[ii]]$bionScore <- bionScore
  msms_sel[[ii]]$yionScore <- yionScore
  msms_sel[[ii]]$ricScore <- ricScore
  msms_sel[[ii]]$bonds <- bonds
  msms_sel[[ii]]$ricScoreRemark <- remark
}

# remove big tables from ic-score script
rm(msms)