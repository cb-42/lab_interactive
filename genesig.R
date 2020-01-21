# Gene expression significance testing script: Christopher Brown - released on 1/21/2020

genesig <- function(infile) {
  # Prompt user to input the significance threshold
  input <- ""
  while(!is.numeric(input)) {
    input <- readline("Enter significance threshold (between 0 and 1; ex: 0.05): ")
    if(grepl("^0?\\.{1}[0-9]+$", input)) {
      sig <- as.numeric(input) # coerce to numeric if only digits (0.0, etc) are present
      break # exit loop
    }
  }

  # Note: column names (groups) will be de-duplicated on load
  dfgen <- read.csv(infile) # assume input is .csv, with levels in columns and genes in rows
  dfgen_t <- t(dfgen) # transpose the dataframe
  df <- setNames(data.frame(rownames(dfgen_t)[-1], dfgen_t[-1,], row.names = NULL, stringsAsFactors = F), c("Group", dfgen_t[1,])) # create dataframe, exclude top row (colnames); add colnames: Group + all gene names
  # create factor; truncate digits (which were created when colnames were de-duplicated on import)
  df$Group <- factor(gsub(gsub(x = colnames(dfgen)[-1], pattern = " ", replacement = "_"), pattern=".\\d", replacement = ""))
  
  # create dataframe containing all gene expression significance tests
  siglist <- vector("list", dim(df)[2]-1)
  for(i in 1:length(siglist)) {
    results <- TukeyHSD(aov(as.numeric(df[,i+1]) ~ df[,"Group"], data = df)) # col needs to be i + 1 to skip group row
    tmpdf <- data.frame(t(results[[1]]), Gene = colnames(df)[i+1], row.names = NULL)[4,] # transpose, keep only "p adj" stat
    ncomps <- sum(1:length(levels(df$Group)) - 1) # number of comparison columns to expect: sum(1:nlevels-1)
    siglist[[i]] <- tmpdf[, c(ncomps+1, 1:ncomps)] # put Gene column first
  }
  resdf <- do.call(rbind, siglist) # reduces list structure to dataframe

  if(ncomps==1) { # Yes if significant
    yndf <- data.frame(ifelse(resdf[,2] < sig, "Yes", "No")) # need to recreate as a dataframe due to 1 column being simplified to vector
  } else {
    yndf <- ifelse(resdf[, 2:(ncomps+1)] < sig, "Yes", "No")
  }
  colnames(yndf) <- paste0(colnames(resdf)[-1], "_signif")
  resdf <- cbind(resdf, yndf)
  
  outfile <- paste0(gsub(x = infile, ".csv", ""), "_stats.csv")
  write.csv(x = resdf, file = outfile, row.names = F) # write results to outfile based on input file
}
