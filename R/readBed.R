# read and validate bed file
readBed <- function(bedFile, bedHeader = FALSE) {
    # read file first
    bed <- data.table::fread(file = bedFile, quote = "", header = bedHeader, 
                             data.table = FALSE, fill = TRUE,
                             stringsAsFactors = FALSE)
    # check number of columns first
    if (ncol(bed) < 3) {
      stop(paste0('Invaild bed file. Please provide vaild bed file that has',
                  ' at least 3 columns.'))
      # position should be integers
    }else if (typeof(bed[, 2]) != 'integer' | typeof(bed[, 3]) != 'integer') {
      # if character only exists in the first 
      stop(paste0('Invaild bed file. 2nd and 3rd columns should be integers.',
                  ' Maybe your bed file has header?'))
    }else{
      # add "chr" if it does not exist in the first column
      if (!(str_detect(bed[, 1][1], "chr"))) {
        bed[, 1] <- paste0("chr", bed[, 1])
      }
    }
    return(bed)
}
