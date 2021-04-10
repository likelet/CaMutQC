selectTrans <- function(maf, CSQdf, INFOdf) {
  # extract, assign and filter CSQ info from maf data frame
  
  for (n in 1:nrow(maf)) {
    CSQ_general <- strsplit(INFOdf$CSQ[n], split = ",")[[1]]
    CSQ_subinfo <- CSQdf[1:length(CSQ_general), ]
    rownames(CSQ_subinfo) <- 1:length(CSQ_general)
    
    for (m in 1:length(CSQ_general)) {
      CSQ_p <- strsplit(CSQ_general[m], split = "\\|")[[1]]
      if(length(CSQ_p) == (ncol(CSQ_subinfo) - 1)) {
        CSQ_p <- c(CSQ_p, "")
      }
      CSQ_subinfo[m, ] <- CSQ_p
    }
    
    # compare and select
    ## construct the mut data frame
    mutDatFrame <- CSQ_subinfo[, c('BIOTYPE', 'Consequence', 'cDNA_position')]
    for (d in 1:nrow(mutDatFrame)) {
      for (e in 1:ncol(mutDatFrame)) {
        if (mutDatFrame[d, e] == '')
          mutDatFrame[d, e] <- 'Missing'
      }
    }
    
    CSQnum <- selectMut(mutDatFrame)
    CSQdf[n, ] <- CSQ_subinfo[CSQnum, ]
  }
  
  return(CSQdf)
}
