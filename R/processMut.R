#' processMut
#' @description Takes union or intersection on multiple MAF data frame, 
#' and return 7 important columns.
#' @param mafList A list of MAF data frames after going through
#' at least one CaMutQC filtration function, and the length of the list <= 3.
#' @param processMethod Methods for processing mutations, including "union"
#' and "intersection". Default: "union".
#' @return A data frame includes mutations after taking union or intersection.
#' @import dplyr
#' @importFrom methods is
#'
#' @export processMut
#' @examples
#' maf_MuSE <- vcfToMAF(system.file("extdata/Multi-caller",
#' "WES_EA_T_1.MuSE.vep.vcf", package="CaMutQC"))
#' maf_MuSE_f <- mutFilterCom(maf_MuSE, report=FALSE, TMB=FALSE,
#' PONfile=system.file("extdata"," PON_test.txt", package="CaMutQC"), 
#' PONformat="txt")
#' maf_VarScan2 <- vcfToMAF(system.file("extdata/Multi-caller",
#' "WES_EA_T_1_varscan_filter_snp.vep.vcf", package="CaMutQC"))
#' maf_VarScan2_f <- mutFilterCom(maf_VarScan2, report=FALSE, TMB=FALSE,
#' PONfile=system.file("extdata","PON_test.txt", package="CaMutQC"), 
#' PONformat="txt")
#' mafs <- list( maf_MuSE_f, maf_VarScan2_f)
#' maf_union <- processMut(mafs, processMethod="union")

processMut <- function(mafList, processMethod = "union") {
    ## check user input, and remove labeled mutations
    filtered_mafs <- list()
    for (i in seq_len(length(mafList))){
        if (!('CaTag' %in% colnames(mafList[[i]])) | !(is(mafList[[i]], "data.frame"))){
            mes <- paste0("'CaTag' column is missing in the ", i, " MAF frame, ",
                          "please make sure all input MAF have gone through at least", 
                          " 1 CaMutQC filtration function.")
            stop(mes)
        }else{
            filtered_mafs[[i]] <- mafList[[i]][which(mafList[[i]]$CaTag == '0'), ]
        }
    }
    ## take union or intersection on filtered mafs
    if ((processMethod == "union") || (processMethod == "intersection")){
        return(processMafs(filtered_mafs, processMethod))
    }else{
        stop("Invalid process method, shoule be 'union' or 'intersection'.")
    }
}

# get the unio or intersection 
processMafs <- function(mafs, method){
    # define base case
    if (length(mafs) == 0){
        stop("No MAF data frame detected in the list!")
    }else if (length(mafs) == 1){
        return(mafs[[1]])
    }else if (length(mafs) == 2){
        if (method == "union"){
            return(dplyr::union(mafs[[1]][, c(1, 4, 5, 6, 7, 9, 10, 11, 13)], 
                                    mafs[[2]][, c(1, 4, 5, 6, 7, 9, 10, 11, 13)]))
        }else{
            return(dplyr::intersect(mafs[[1]][, c(1, 4, 5, 6, 7, 9, 10, 11, 13)], 
                                    mafs[[2]][, c(1, 4, 5, 6, 7, 9, 10, 11, 13)]))
        }
    }else if (length(mafs) == 3){
        if (method == "union"){
            temp_maf <- dplyr::union(mafs[[1]][, c(1, 4, 5, 6, 7, 9, 10, 11, 13)], 
                                    mafs[[2]][, c(1, 4, 5, 6, 7, 9, 10, 11, 13)])
            return(dplyr::union(temp_maf, mafs[[3]][, c(1, 4, 5, 6, 7, 9, 10, 11, 13)]))
        }else{
            temp_maf <- dplyr::intersect(mafs[[1]][, c(1, 4, 5, 6, 7, 9, 10, 11, 13)], 
                                    mafs[[2]][, c(1, 4, 5, 6, 7, 9, 10, 11, 13)])
            return(dplyr::intersect(temp_maf, mafs[[3]][, c(1, 4, 5, 6, 7, 9, 10, 11, 13)]))
        }
    }else {
        stop("Too many mafs, please provide at most 3 maf data frames.")
    }
}
