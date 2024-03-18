selectTrans <- function(csqDat, infoDat) {
    # extract, assign and filter CSQ info from infoDat data frame
    for (n in seq_len(nrow(csqDat))) {
      csqGeneral <- strsplit(infoDat$CSQ[n], split=",")[[1]]
      # split the csq info 
      csqSubInfo <- data.frame(matrix(unlist(lapply(csqGeneral, batchCSQ, csqDat)), 
                              nrow=length(csqGeneral), byrow=TRUE), 
                              stringsAsFactors=FALSE)
      rownames(csqSubInfo) <- seq_len(length(csqGeneral))
      colnames(csqSubInfo) <- colnames(csqDat)
      # compare and select
      ## construct the mut data frame
      mutDatFrame <- csqSubInfo[, c('BIOTYPE', 'Consequence', 'cDNA_position')]
      mutDatFrame[mutDatFrame == "" | grepl('\\?', mutDatFrame)] <- "Missing"
      CSQnum <- selectMut(mutDatFrame)
      csqDat[n, ] <- csqSubInfo[CSQnum, ]
    }
    return(csqDat)
}

# helper function for splitting CSQ info
batchCSQ <- function(csq, csqDat) {
    csqSingle <- strsplit(csq, split="\\|")[[1]]
    if (length(csqSingle) == (ncol(csqDat) - 1)) {
        csqSingle <- c(csqSingle, "")
    }
    return(csqSingle)
}

