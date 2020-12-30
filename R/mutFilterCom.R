#' mutFilterCom
#' @description Apply common filter strategies on a MAF data frame.
#'
#' @param maf An MAF data frame.
#' @param tumorSampleName Label of the tumor sample, should be one of the 
#' column names of maf. If it is set as 'Extracted', tumorSampleName would be
#' extracted automatically from the maf data frame. Default: 'Extracted'.
#' @param normalSampleName Label of the normal sample, should be one of the 
#' column names of maf. If it is set as 'Extracted', normalSampleName would be
#' extracted automatically from the maf data frame. Default: 'Extracted'.
#' @param filter Whether to directly return a filtered MAF data frame. If FALSE,
#' a simulation filtering process will be run, 
#' and a filter report may be generated.
#' @param show.report Whether to return the filtering report. Default: TRUE
#' 
#' @return An MAF data frame after SNP filtering

mutFilterCom <- function(maf, tumorSampleName = 'Extracted', 
                         normalSampleName = 'Extracted', filter = FALSE, 
                         show.report = TRUE) {
  
  
}
