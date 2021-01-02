filterReport <- function(oriMAF, filtMAF) {
  
  # check status
  if (nrow(oriMAF) == 0){
    stop('No information found in oriMAF.')
  }else if(nrow(filtMAF) == 0){
    return('No variants left after filter.')
  }
  
  # extract info contained in the report
  ## number of variants
  nVarOri <- nrow(oriMAF)
  nVarFilt <- nrow(filtMAF)
  
  ## number of genes
  nGeneOri <- length(unique(oriMAF$Hugo_Symbol))
  nGeneFilt <- length(unique(filtMAF$Hugo_Symbol))
  
  ## location in genome (spread in chromosome)
  locOri <- paste(names(table(oriMAF$Chromosome)), 
                  as.character(table(oriMAF$Chromosome)), collapse = "  ", 
                  sep = ":")
  locFilt <- paste(names(table(filtMAF$Chromosome)), 
                   as.character(table(filtMAF$Chromosome)), collapse = "  ", 
                   sep = ":")
  
  
  ## variant classification
  vclaOri <- paste(names(table(oriMAF$Variant_Classification)), 
                  as.character(table(oriMAF$Variant_Classification)), 
                  collapse = "  ", sep = ":")
  vclaFilt <- paste(names(table(filtMAF$Variant_Classification)), 
                   as.character(table(filtMAF$Variant_Classification)), 
                   collapse = "  ", sep = ":")
  
  ## variant types
  vtypeOri <- paste(names(table(oriMAF$Variant_Type)), 
                  as.character(table(oriMAF$Variant_Type)), 
                  collapse = "  ", sep = ":")
  vtypeFilt <- paste(names(table(filtMAF$Variant_Type)), 
                   as.character(table(filtMAF$Variant_Type)), 
                   collapse = "  ", sep = ":")
  
  ## average VAF
  VAFori <- mean(oriMAF$VAF)
  VAFFilt <- mean(filtMAF$VAF)
  
  
  # export report
  report <- as.data.frame(matrix(ncol = 6, nrow = 2))
  colnames(report) <- c('#Variant', '#Gene', 'Location', 
                              'Variant_Classification', 'Variant_Type', 
                              'VAFMean')
  rownames(report) <- c('Before', 'After')
  report[1, 1] <- nVarOri
  report[2, 1] <- nVarFilt
  report[1, 2] <- nGeneOri
  report[2, 2] <- nGeneFilt
  report[1, 3] <- locOri
  report[2, 3] <- locFilt
  report[1, 4] <- vclaOri
  report[2, 4] <- vclaFilt
  report[1, 5] <- vtypeOri
  report[2, 5] <- vtypeFilt
  report[1, 6] <- VAFori
  report[2, 6] <- VAFFilt
  return(report)
}


