#' mutFilterCom
#' @description Apply common filter strategies on a MAF data frame.
#'
#' @param maf An MAF data frame.
#' @param tumorSampleName Label of the tumor sample, should be one of the 
#' column names of MAF. If it is set as 'Extracted', tumorSampleName would be
#' extracted automatically from the MAF data frame. Default: 'Extracted'.
#' @param normalSampleName Label of the normal sample, should be one of the 
#' column names of MAF. If it is set as 'Extracted', normalSampleName would be
#' extracted automatically from the MAF data frame. Default: 'Extracted'.
#' @param mut.filter Whether to directly return a filtered MAF data frame. 
#' If FALSE, a simulation filtering process will be run,  and only a filter 
#' report will be generated. If TRUE, a new MAF data frame and a filter report 
#' will be generated. Default: TRUE
#' @param report.dir Path to the output report file. Default: './'
#' @import rmarkdown
#' 
#' @return An MAF data frame after SNP filtering

mutFilterCom <- function(maf, tumorSampleName = 'Extracted', 
                         normalSampleName = 'Extracted', mut.filter = TRUE,
                         report.dir = './') {
  
  # process tumorSampleName and normalSampleName
  if (tumorSampleName == 'Extracted'){
    tumorSample <- unique(maf$Tumor_Sample_Barcode)
  } else if(!(tumorSampleName %in% colnames(maf))) {
    stop('Invaild tumorSampleName.')
  }
  
  if (normalSampleName == 'Extracted'){
    normalSample <- unique(maf$Matched_Norm_Sample_Barcode)
  }else if(!(normalSampleName %in% colnames(maf))) {
    stop('Invaild normalSampleName.')
  }
 
  
  # start filtering 
  ## sequencing quality
  cat('Variants filtering for sequencing quality is in process.')
  maf_q <- mutFilterQual(maf, tumorSampleName = tumorSample, 
                         normalSampleName = normalSample)
  if (nrow(maf_q) == 0) {
    message('No variation left after filtering variants in low sequencing quality.')
  }
  
  ## PON
  cat('\nVariants filtering for PON is in process.\n')
  maf_p <- mutFilterPON(maf_q)
  if (nrow(maf_p) == 0) {
    message('No variation left after filtering variants in PON.')
  }
  
  ## SNP
  cat('\nVariants filtering of SNP variants is in process.')
  maf_snp <- mutFilterSNP(maf_p)
  if (nrow(maf_snp) == 0) {
    message('No variation left after filtering SNP variants.')
  }
  
  ## variant types
  cat('\nVariants filtering of variant types.\n')
  maf_type <- mutFilterType(maf_snp)
  if (nrow(maf_type) == 0) {
    message('No variation left after filtering variants of certain types.')
  }
  
  ## Strand Bias
  cat('\nVariants filtering of strand bias is in process.\n')
  #cat(tumorSample %in% colnames(maf_type))
  maf_sb <- mutFilterSB(maf_type, tumorSampleName = tumorSample)
  if (nrow(maf_sb) == 0) {
    message('No variation left after filtering SB variants.')
  }
  
  ## normal depth
  cat('\nVariants filtering of normal depth is in process.')
  maf_n <- mutFilterNormalDP(maf_sb, normalSampleName = normalSample)
  if (nrow(maf_n) == 0) {
    message('No variation left after filtering variants through normal depth.')
  }
  
  ## COSMIC
  cat('\nVariants filtering of non-COSMIC variants is in process.\n')
  maf_c <- mutFilterCOSMIC(maf_n)
  if (nrow(maf_c) == 0) {
    message('No variation left after filtering variants not in COSMIC.')
  }
  
  # report generation
  if (mut.filter) {
    rmarkdown::render('./report/report.Rmd', paste0(report.dir, 
                                            'filterReport.html'))
    return(maf_c)
  }else{
    rmarkdown::render('./report/report.Rmd', paste0(report.dir, 
                                            'filterReport.html'))
  }
}
