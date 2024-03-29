readINFO <- function(vcfHeader, vcfMain) {
    Info <- strsplit(vcfHeader[,1][grep("INFO=<ID=", vcfHeader[,1])],
                     split="##INFO=<ID=")
    for (i in seq_len(length(Info))) {
        Info[[i]] <- Info[[i]][2]
    }
    Info <- as.character(Info)
    # construct the info frame
    Info_frame <- data.frame(matrix(ncol=4, nrow=length(Info)))
    colnames(Info_frame) <- c('ID', 'Number', 'Type', 'Description')
    for (j in seq_len(length(Info))){
      Info_frame[j, 1] <- strsplit(Info[j], split=",", fixed=TRUE)[[1]][1]
      Info_frame[j, 2] <- strsplit(strsplit(Info[j], split=",", 
                                            fixed=TRUE)[[1]][2], split="=")[[1]][2]
      Info_frame[j, 3] <- strsplit(strsplit(Info[j], split=",", 
                                            fixed=TRUE)[[1]][3], split="=")[[1]][2]
      Info_frame[j, 4] <- strsplit(strsplit(Info[j], split=",", 
                                            fixed=TRUE)[[1]][4], split="=")[[1]][2]
      Info_frame[j, 4] <- str_sub(Info_frame[j, 4], start=2,
                                  end=nchar(Info_frame[j, 4])-2)
    }
    # construct Info data frame
    Info_m <- as.data.frame(matrix(ncol=nrow(Info_frame), 
                                   nrow=nrow(vcfMain)))
    colnames(Info_m) <- Info_frame$ID
    ## extract data from Info column and fill in blanks in Info_m frame
    for(r in seq_len(nrow(Info_m))) {
      Infos <- strsplit(vcfMain[r, 8], split=";", fixed=TRUE)[[1]]
      for (t in seq_len(length(Infos))) {
        ID <- str_split(Infos[t], pattern='(?<=[:alpha:])\\=')[[1]][1]
        inform <- str_split(Infos[t], pattern='(?<=[:alpha:])\\=')[[1]][2]
        if(ID %in% colnames(Info_m)) {
            Info_m[r, ID] <- inform
        }
      }
    }
    return(list(Info_frame, Info_m))
}
