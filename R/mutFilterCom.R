#' mutFilterCom
#' @description Apply common filter strategies on a MAF data frame.
#'
#' @param maf An MAF data frame.
#' @param tumorDP Threshold of tumor total depth. Default: 20
#' @param normalDP Threshold of normal total depth. Default: 10
#' @param tumorAD Threshold of tumor alternative allele depth. Default: 10
#' @param VAF Threshold of VAF value. Default: 0.05
#' @param VAFratio Threshold of VAF ratio (tVAF/nVAF). Default: 5
#' @param tumorSampleName Label of the tumor sample, should be one of the 
#' column names of maf. If it is set as 'Extracted', tumorSampleName would be
#' extracted automatically from the maf data frame. Default: 'Extracted'.
#' @param SBmethod Method will be used to detect strand bias, 
#' including 'SOR' and 'Fisher'. Default: 'SOR'. SOR: StrandOddsRatio 
#' (https://gatk.broadinstitute.org/hc/en-us/articles/360041849111-
#' StrandOddsRatio)
#' @param SBscore Cutoff strand bias score used to filter variants.
#' Default: 3
#' @param maxIndelLen Maximum length of indel accepted to be included. 
#' Default: 50
#' @param minInterval Maximum length of interval between an SNV and an indel 
#' accepted to be included. Default: 10
#' @param tagFILTER Variants with spcific tag in the FILTER column will be kept,
#' Default: 'PASS'
#' @param dbVAF Threshold of VAF of certain population for variants
#'  in database. Default: 0.01
#' @param ExAC Whether to filter variants listed in ExAC with VAF higher than
#' cutoff(set in VAF parameter). Default: TRUE.
#' @param Genomesprojects1000 Whether to filter variants listed in 
#' Genomesprojects1000 with VAF higher than cutoff(set in VAF parameter). 
#' Default: TRUE.
#' @param ESP6500 Whether to filter variants listed in ESP6500 with VAF higher 
#' than cutoff(set in VAF parameter). Default: TRUE.
#' @param gnomAD Whether to filter variants listed in gnomAD with VAF higher 
#' than cutoff(set in VAF parameter). Default: TRUE.
#' @param dbSNP Whether to filter variants listed in dbSNP. Default: TRUE.
#' @param COSMIConly Whether to only keep variants in COSMIC. Default: FALSE.
#' @param keepType A group of variant classifications will be kept, 
#' including 'exonic' and 'nonsynonymous'. Default: 'exonic'. 
#' @param bedFile A file in bed format that contains region information.
#' Default: NULL
#' @param bedFilter Whether to filter the information in bed file or not, which 
#' only leaves segments in Chr1-Ch22, ChrX and ChrY. Default: TRUE
#' @param mutFilter Whether to directly return a filtered MAF data frame. 
#' If FALSE, a simulation filtration process will be run, and the original MAF 
#' data frame with tags in CaTag column, and  a filter report will be returned. 
#' If TRUE, a filtered MAF data frame and a filter report will be generated. 
#' Default: FALSE
#' @param selectCols Columns will be contained in the filtered data frame.
#' By default (TRUE), the first 13 columns and 'Tumor_Sample_Barcode' column.
#' Or a vector contains column names will be kept.
#' @param report Whether to generate report automatically. Default: TRUE
#' @param reportFile File name of the report. Default: 'FilterReport.html'
#' @param reportDir Path to the output report file. Default: './'
#' @param TMB Whether to calculate TMB. Default: TRUE
#' 
#' @return An MAF data frame after common strategy filtration
#' @return A filter report in HTML format
#' 
#' @export mutFilterCom

mutFilterCom <- function(maf, tumorDP = 20, normalDP = 10, tumorAD = 10, 
                         VAF = 0.05, VAFratio = 5, tumorSampleName = 'Extracted', 
                         SBmethod = 'SOR', SBscore = 3, maxIndelLen = 50, 
                         minInterval = 10, tagFILTER = 'PASS', dbVAF = 0.01, 
                         ExAC = TRUE, Genomesprojects1000 = TRUE, ESP6500 = TRUE, 
                         gnomAD = TRUE, dbSNP = TRUE, COSMIConly = TRUE, 
                         keepType = 'exonic', bedFile = NULL, bedFilter = TRUE, 
                         mutFilter = FALSE, selectCols = TRUE, report = TRUE, 
                         reportFile = 'FilterReport.html', reportDir = './', 
                         TMB = TRUE) {
  
  # process tumorSampleName
  if (tumorSampleName == 'Extracted'){
    tumorSample <- unique(maf$Tumor_Sample_Barcode)
  } else if(!(tumorSampleName %in% colnames(maf))) {
    stop('Invaild tumorSampleName.')
  }
  
  # run mutFilterTech
  mafFilteredT <- mutFilterTech(maf, tumorDP = tumorDP, normalDP = normalDP, 
                               tumorAD = tumorAD, VAF = VAF, VAFratio = VAFratio, 
                               tumorSampleName = tumorSampleName, 
                               SBmethod = SBmethod, SBscore = SBscore, 
                               maxIndelLen = maxIndelLen, 
                               minInterval = minInterval, tagFILTER = tagFILTER)
  
  # filter first for report usage
  mafFilteredTs <- mafFilteredT[mafFilteredT$CaTag == '0', ]
  
  # run mutSelection
  mafFilteredS <- mutSelection(mafFilteredT, dbVAF = dbVAF, ExAC = ExAC, 
                              Genomesprojects1000 = Genomesprojects1000, 
                              ESP6500 = ESP6500, gnomAD = gnomAD, dbSNP = dbSNP,
                              COSMIConly = COSMIConly, keepType = keepType,
                              bedFile = bedFile, bedFilter = bedFilter)
  
  # filter first for report usage
  mafFilteredS2 <- suppressMessages(
    mutSelection(mafFilteredTs, dbVAF = dbVAF, ExAC = ExAC, 
                 Genomesprojects1000 = Genomesprojects1000, dbSNP = dbSNP,
                 ESP6500 = ESP6500, gnomAD = gnomAD, COSMIConly = COSMIConly, 
                 keepType = keepType, bedFile = bedFile, bedFilter = bedFilter))
  
  mafFilteredF <- mafFilteredS2[mafFilteredS2$CaTag == '0', ]
  # print(nrow(mafFilteredF))
  if (nrow(mafFilteredF) == 0){
    stop('No variants left after filtration.')
  }
  
  if (TMB){
    # check bed file
    if (is.null(bedFile)){
      stop(paste0('A bed file is missing, which is required for TMB calculation.',
          ' If you don\'t want to calculate TMB, please set TMB to FALSE.'))
    }else{
      bed <- read.table(bedFile)
      bedLen <- as.character(round(sum(bed$V3 - bed$V2)/1000000, 2))
      TMBvalue <- calTMB(mafFilteredF, bedFile = bedFile)
    }
  }
  
  # report generation
  if (report){
    rmarkdown::render('./report/FilterReport.Rmd', output_file = reportFile,
                      output_dir = reportDir)
  }
  
  if (mutFilter) {
    if (selectCols){
      if (isTRUE(selectCols)){
        return(mafFilteredF[, c(1:12, 16)])
      }else{
        if (all(selectCols %in% colnames(mafFilteredF))){
          return(mafFilteredF[, selectCols])
        }else{
          stop('Not all selected columns can be found in MAF columns. ')
        }
      }
    }
  }else{
    return(mafFilteredS)
  }
}
