#' mutFilterNormalDP
#' @description Filter dbsnp/non-dbsnp variants based on their normal depth.
#' Variants in dbSNP database should have normal depth >= 19, while non-dbSNP
#' variants should have normal depth >= 8 to avoid being filtered.
#' @param maf An MAF data frame, generated by \code{\link{vcfToMAF}} function.
#' @param dbsnpCutoff Cutoff of normal depth for dbSNP variants. Default: 19.
#' @param nonCutoff Cutoff of normal depth for non-dbSNP variants. Default: 8.
#' @param verbose Whether to generate message/notification during the 
#' filtration process. Default: TRUE.
#' @importFrom methods is
#' @import dplyr
#' 
#' @return An MAF data frame where some variants
#' has N tag in CaTag column for Normal depth filtration.
#' @export mutFilterNormalDP
#' @examples
#' maf <- vcfToMAF(system.file("extdata",
#' "WES_EA_T_1_mutect2.vep.vcf", package="CaMutQC"))
#' mafF <- mutFilterNormalDP(maf)


mutFilterNormalDP <- function(maf, dbsnpCutoff = 19, nonCutoff = 8, 
                              verbose = TRUE) {
    # check user input
    if (!(is(maf, "data.frame"))) {
      stop("maf input should be a data frame, did you get it from vcfToMAF function?")
    }
    # filter based on normal DP
    if (!('Existing_variation' %in% colnames(maf))) {
      discard <- c()
      if (verbose) {
          mes <- paste0('No annotation from existing database detected, ', 
                        "mutFilterNormalDP function will not run.")
          message(mes)
      }
    }
    # build a vector to store discarded variants
    discard <- maf %>% 
      # convert to tibble and name the rownames column
      as_tibble(rownames="rowname") %>% 
      filter((grepl('rs', Existing_variation) & n_depth < dbsnpCutoff) | 
             (!grepl('rs', Existing_variation) & n_depth < nonCutoff)) %>%
    pull(rowname)
    # add tag
    discard <- as.vector(na.omit(discard))
    maf[discard, 'CaTag'] <- paste0(maf[discard, 'CaTag'], 'N')
    return(maf)
}

