# read and validate bed file / bed data frame or bed rds
readBed <- function(bedVar, bedHeader = FALSE) {
    if (is(bedVar, "data.frame")) {
        bed <- bedVar
    # read it if it is a bed file
    }else if (substr(bedVar, nchar(bedVar)-2, nchar(bedVar)) == "bed"){
        # read file first
        bed <- data.table::fread(file=bedVar, quote="", header=bedHeader, 
                               data.table=FALSE, fill=TRUE,
                               stringsAsFactors=FALSE)
    }else if (substr(bedVar, nchar(bedVar)-2, nchar(bedVar)) == "rds"){
        # load rds file
        load(bedVar)
    }else{
        # check whether it is null or invalid
        stop('Please provide a bed file, a dataframe or a rds file.')
    }
    # check number of columns first
    if (ncol(bed) < 3) {
        mes <- paste0('Invaild bed file. Please provide vaild bed file that has',
                      ' at least 3 columns.')
        stop(mes)
      # position should be integers
    }else if (typeof(bed[, 2]) != 'integer' | typeof(bed[, 3]) != 'integer') {
        # if character only exists in the first 
        mes <- paste0('Invaild bed file. 2nd and 3rd columns should be integers.',
                        ' Maybe your bed file has headers?')
        stop(mes)
    }else{
        # add "chr" if it does not exist in the first column
        if (any(!str_detect(bed[, 1], "chr"))) {
            noChrRow <- !str_detect(bed[, 1], "chr")
            bed[noChrRow, 1] <- paste0("chr", bed[noChrRow, 1])
        }
    }
    return(bed)
}
