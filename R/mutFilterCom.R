#' mutFilterCom
#' @description Apply common filter strategies on a MAF data frame.
#'
#' @param maf An MAF data frame.
#' @param tumorSampleName Label of the tumor sample, should be one of the 
#' column names of MAF. If it is set as 'Extracted', tumorSampleName would be
#' extracted automatically from the MAF data frame. Default: 'Extracted'.
#' @param mut.filter Whether to directly return a filtered MAF data frame. 
#' If FALSE, a simulation filtration process will be run,  and only a filter 
#' report will be generated. If TRUE, a new MAF data frame and a filter report 
#' will be generated. Default: TRUE
#' @param filterParam Value in FILTER field used to filter variants.
#' Default: 'PASS'
#' @param report.dir Path to the output report file. Default: './'
#' @param TMB Whether to calculate TMB. Default: FALSE
#' @param bedFile A bed File used to calculate TMB. Default: NULL
#' @import rmarkdown, ggplot2
#' 
#' @return An MAF data frame after common strategy filtration
#' @return A filter report in HTML format
#' 
#' @export mutFilterCom

mutFilterCom <- function(maf, tumorSampleName = 'Extracted', mut.filter = TRUE, 
                         filterParam = 'PASS', report.dir = './', TMB = FALSE,
                         bedFile = NULL) {
  
  # process tumorSampleName
  if (tumorSampleName == 'Extracted'){
    tumorSample <- unique(maf$Tumor_Sample_Barcode)
  } else if(!(tumorSampleName %in% colnames(maf))) {
    stop('Invaild tumorSampleName.')
  }
  
  # start filtration 
  ## FILTER field
  cat('Variants filtration based on FILTER field is in process.')
  maf_filtered <- maf[maf$FILTER ==  filterParam, ]
  if (nrow(maf_filtered) == 0) {
    message('No variation left after FILTER field filtration.')
  }else{
    
    ## sequencing quality
    cat('\nVariants filtration for sequencing quality is in process.')
    maf_filtered <- mutFilterQual(maf_filtered)
    if (nrow(maf_filtered) == 0) {
      message('No variation left after filtration variants in low 
              sequencing quality.')
    }else{
      
      ## PON
      cat('\nVariants filtration for PON is in process.\n')
      maf_filtered <- mutFilterPON(maf_filtered)
      if (nrow(maf_filtered) == 0) {
        message('No variation left after filtration variants in PON.')
      }else{
        
        ## SNP
        cat('\nVariants filtration of SNP variants is in process.')
        maf_filtered <- mutFilterSNP(maf_filtered)
        if (nrow(maf_filtered) == 0) {
          message('No variation left after filtration SNP variants.')
        }else{
          ## variant types
          cat('\nVariants filtration of variant types.\n')
          maf_filtered <- mutFilterType(maf_filtered)
          if (nrow(maf_filtered) == 0) {
            message('No variation left after filtration variants of certain types.')
          }else{
            ## Strand Bias
            cat('\nVariants filtration of strand bias is in process.\n')
            #cat(tumorSample %in% colnames(maf_type))
            maf_filtered <- mutFilterSB(maf_filtered, 
                                        tumorSampleName = tumorSample)
            if (nrow(maf_filtered) == 0) {
              message('No variation left after filtration SB variants.')
            }else{
              ## normal depth
              cat('\nVariants filtration of normal depth is in process.')
              maf_filtered <- mutFilterNormalDP(maf_filtered)
              if (nrow(maf_filtered) == 0) {
                message('No variation left after filtration variants through normal depth.')
              }else{
                ## COSMIC
                cat('\nVariants filtration of non-COSMIC variants is in process.\n')
                maf_filtered <- mutFilterCOSMIC(maf_filtered)
                if (nrow(maf_filtered) == 0) {
                  message('No variation left after filtration variants not 
                          in COSMIC.')
                }else{
                  ## mutFilterAdj
                  cat('\nVariants filtration of Adjacent indel is in process.\n')
                  maf_filtered <- mutFilterAdj(maf_filtered)
                }
              }
            }
          }
        }
      }
    }
  }
  
  if (nrow(maf_filtered) == 0){
    stop('No variant left after filtration.')
  }
  
  if (TMB){
    # check bed file
    if (is.null(bedFile)){
      stop('A bed file is missing, which is required for TMB calculation.
           If you don\'t want to calculate TMB, please set TMB as FALSE.')
    }else{
      bed <- read.table(bedFile)
      bedLen <- as.character(round(sum(bed$V3 - bed$V2)/1000000, 2))
      TMBvalue <- calTMB(maf_filtered, bedFile = bedFile)
    }
  }
  
  # report generation
  rmarkdown::render('./report/FilterReport.Rmd', 
                    output_file = paste0(report.dir, 'filterReport.html'))
  if (mut.filter) {
    return(maf_filtered)
  }
}
