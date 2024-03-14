selectTrans <- function(csqDat, infoDat) {
    # extract, assign and filter CSQ info from infoDat data frame
    for (n in seq_len(nrow(csqDat))) {
      csqGeneral <- strsplit(infoDat$CSQ[n], split=",")[[1]]
      csqSubInfo <- csqDat[seq_len(length(csqGeneral)), ]
      rownames(csqSubInfo) <- seq_len(length(csqGeneral))
  
      for (m in seq_len(length(csqGeneral))) {
          csqSingle <- strsplit(csqGeneral[m], split="\\|")[[1]]
          if(length(csqSingle) == (ncol(csqSubInfo) - 1)) {
            csqSingle <- c(csqSingle, "")
          }
          csqSubInfo[m, ] <- csqSingle
      }
      # compare and select
      ## construct the mut data frame
      mutDatFrame <- csqSubInfo[, c('BIOTYPE', 'Consequence', 'cDNA_position')]
      for (d in seq_len(nrow(mutDatFrame))) {
          for (e in seq_len(ncol(mutDatFrame))) {
            if ((mutDatFrame[d, e] == '') | (grepl('\\?', mutDatFrame[d, e])))
              mutDatFrame[d, e] <- 'Missing'
          }
      }
      CSQnum <- selectMut(mutDatFrame)
      csqDat[n, ] <- csqSubInfo[CSQnum, ]
    }
    return(csqDat)
}
