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
        # first replace empty elements into "missing"
        mutDatFrame[mutDatFrame == ""] <- "Missing"
        # then replace elements with "?" into missing
        indices <- findQuestionMarks(mutDatFrame)
        mutDatFrame[indices$BIOTYPE, 1] <- "Missing"
        mutDatFrame[indices$Consequence, 2] <- "Missing"
        mutDatFrame[indices$cDNA_position, 3] <- "Missing"
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

# helper function to help find elements with "?" in a data frame
findQuestionMarks <- function(dataframe) {
    indices <- sapply(dataframe, function(col) grep("\\?", col))
    return(indices)
}

