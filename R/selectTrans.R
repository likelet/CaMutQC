selectTrans <- function(maf, csqDat, infoDat) {
    # extract, assign and filter CSQ info from maf data frame
    for (n in seq_len(nrow(maf))) {
      CSQ_general <- strsplit(infoDat$CSQ[n], split=",")[[1]]
      CSQ_subinfo <- csqDat[seq_len(length(CSQ_general)), ]
      rownames(CSQ_subinfo) <- seq_len(length(CSQ_general))
  
      for (m in seq_len(length(CSQ_general))) {
          CSQ_p <- strsplit(CSQ_general[m], split="\\|")[[1]]
          if(length(CSQ_p) == (ncol(CSQ_subinfo) - 1)) {
            CSQ_p <- c(CSQ_p, "")
          }
          CSQ_subinfo[m, ] <- CSQ_p
      }
      # compare and select
      ## construct the mut data frame
      mutDatFrame <- CSQ_subinfo[, c('BIOTYPE', 'Consequence', 'cDNA_position')]
      for (d in seq_len(nrow(mutDatFrame))) {
          for (e in seq_len(ncol(mutDatFrame))) {
            if ((mutDatFrame[d, e] == '') | (grepl('\\?', mutDatFrame[d, e])))
              mutDatFrame[d, e] <- 'Missing'
          }
      }
      CSQnum <- selectMut(mutDatFrame)
      csqDat[n, ] <- CSQ_subinfo[CSQnum, ]
    }
    return(csqDat)
}
