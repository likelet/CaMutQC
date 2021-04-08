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
#' @param filterParam Value in FILTER field used to filter variants.
#' Default: 'PASS'
#' @param report.dir Path to the output report file. Default: './'
#' @param TMB Whether to calculate TMB. Default: FALSE
#' @param bedFile A bed File used to calculate TMB. Default: NULL
#' @import rmarkdown, ggplot2
#' 
#' @return An MAF data frame after common strategy filtering
#' @return A filter report in HTML format
#' 
#' @export mutFilterCom

mutFilterCom <- function(maf, tumorSampleName = 'Extracted', 
                         normalSampleName = 'Extracted', mut.filter = TRUE, 
                         filterParam = 'PASS', report.dir = './', TMB = FALSE,
                         bedFile = NULL) {
  
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
  ## FILTER field
  cat('Variants filtering based on FILTER field is in process.')
  maf_filtered <- maf[maf$FILTER ==  filterParam, ]
  if (nrow(maf_filtered) == 0) {
    message('No variation left after FILTER field filtering.')
  }else{
    
    ## sequencing quality
    cat('\nVariants filtering for sequencing quality is in process.')
    maf_filtered <- mutFilterQual(maf_filtered, tumorSampleName = tumorSample, 
                           normalSampleName = normalSample)
    if (nrow(maf_filtered) == 0) {
      message('No variation left after filtering variants in low 
              sequencing quality.')
    }else{
      
      ## PON
      cat('\nVariants filtering for PON is in process.\n')
      maf_filtered <- mutFilterPON(maf_filtered)
      if (nrow(maf_filtered) == 0) {
        message('No variation left after filtering variants in PON.')
      }else{
        
        ## SNP
        cat('\nVariants filtering of SNP variants is in process.')
        maf_filtered <- mutFilterSNP(maf_filtered)
        if (nrow(maf_filtered) == 0) {
          message('No variation left after filtering SNP variants.')
        }else{
          ## variant types
          cat('\nVariants filtering of variant types.\n')
          maf_filtered <- mutFilterType(maf_filtered)
          if (nrow(maf_filtered) == 0) {
            message('No variation left after filtering variants of certain types.')
          }else{
            ## Strand Bias
            cat('\nVariants filtering of strand bias is in process.\n')
            #cat(tumorSample %in% colnames(maf_type))
            maf_filtered <- mutFilterSB(maf_filtered, 
                                        tumorSampleName = tumorSample)
            if (nrow(maf_filtered) == 0) {
              message('No variation left after filtering SB variants.')
            }else{
              ## normal depth
              cat('\nVariants filtering of normal depth is in process.')
              maf_filtered <- mutFilterNormalDP(maf_filtered, 
                                                normalSampleName = normalSample)
              if (nrow(maf_filtered) == 0) {
                message('No variation left after filtering variants through normal depth.')
              }else{
                ## COSMIC
                cat('\nVariants filtering of non-COSMIC variants is in process.\n')
                maf_filtered <- mutFilterCOSMIC(maf_filtered)
                if (nrow(maf_filtered) == 0) {
                  message('No variation left after filtering variants not 
                          in COSMIC.')
                }
              }
            }
          }
        }
      }
    }
  }
  
  if (TMB){
    # check bed file
    if (is.null(bedFile)){
      stop('A bed file is missing, which is required for TMB calculation.
           If you don\'t want to calculate TMB, please set TMB as FALSE.')
    }else{
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
