#' mutFilterCan
#' @description Apply common filter strategies on a MAF data frame for different
#' cancer types.
#'
#' @param maf An MAF data frame.
#' @param cancerType Type of cancer sample whose params needed to be referred to.
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
#' @return An MAF data frame after common strategy filtration for a cancer type.
#' @return A filter report in HTML format
#' 
#' @export mutFilterCan

mutFilterCan <- function(maf, cancerType, tumorDP = 20, normalDP = 10, 
                         tumorAD = 10, VAF = 0.05, VAFratio = 5, 
                         tumorSampleName = 'Extracted', SBmethod = 'SOR', 
                         SBscore = 3, maxIndelLen = 50, minInterval = 10, 
                         tagFILTER = 'PASS', dbVAF = 0.01, ExAC = TRUE, 
                         Genomesprojects1000 = TRUE, ESP6500 = TRUE, 
                         gnomAD = TRUE, dbSNP = TRUE, COSMIConly = TRUE, 
                         keepType = 'exonic', bedFile = NULL, bedFilter = TRUE, 
                         mutFilter = FALSE, selectCols = TRUE, report = TRUE, 
                         reportFile = 'FilterReport.html', reportDir = './', 
                         TMB = TRUE) {
  
  # turn on the switch for report 
  withType <- TRUE
  # BLCA 
  if (cancerType == 'BLCA'){
    mafFiltered <- mutfilterCom(maf, SBmethod = 'Fisher', SBscore = 20, 
                                minInterval = 30, dbVAF = 0, tumorDP = 10,
                                Genomesprojects1000 = TRUE, ExAC = TRUE, 
                                normalDP = 10, tumorAD = 5, VAF = 0,
                                VAFratio = 0, maxIndelLen = Inf, keepType = 'ALL',
                                tagFILTER = NULL, Genomesprojects1000 = TRUE, 
                                ESP6500 = FALSE, gnomAD = FALSE, dbSNP = FALSE)
  # BRCA
  }else if(cancerType == 'BRCA'){
    mafFiltered <- mutFilterCom(maf, tumorAD = 5, VAF = 0.1, dbVAF = 0,
                                dbSNP = TRUE, Genomesprojects1000 = TRUE,
                                ESP6500 = TRUE, COSMIConly = TRUE, SBscore = 0,
                                keepType = 'ALL', tumorDP = 6, normalDP = 6, 
                                VAFratio = 0, maxIndelLen = Inf, minInterval = 0, 
                                ExAC = FALSE, gnomAD = FALSE, tagFILTER = NULL)
  # COADREAD
  }else if(cancerType == 'COADREAD'){
    mafFiltered <- mutFilterCom(maf, tumorDP = 5, VAF = 0.2, dbVAF = 0,
                                dbSNP = TRUE, Genomesprojects1000 = TRUE,
                                normalDP = 0, tumorAD = 0, VAFratio = 0, 
                                maxIndelLen = Inf, minInterval = 0, SBscore = 0,
                                tagFILTER = NULL, ExAC = FALSE, gnomAD = FALSE,
                                ESP6500 = FALSE, COSMIConly = FALSE)
  # UCEC
  }else if(cancerType == 'UCEC'){
    mafFiltered <- mutFilterCom(maf, tumorAD = 0, tumorDP = 0, normalDP = 0,
                                dbVAF = 0, dbSNP = TRUE, COSMIConly = FALSE,
                                Genomesprojects1000 = TRUE, VAF = 0, SBscore = 0,
                                VAFratio = 0, maxIndelLen = Inf, minInterval = 0, 
                                tagFILTER = NULL, ExAC = FALSE, ESP6500 = FALSE, 
                                gnomAD = FALSE, dbSNP = FALSE)
  
  # UCS
  }else if(cancerType == 'UCS'){
    mafFiltered <- mutFilterCom(maf, tumorAD = 5, tumorDP = 12, normalDP = 5,
                                dbVAF = 0, dbSNP = FALSE, COSMIConly = TRUE,
                                Genomesprojects1000 = FALSE, VAF = 0, 
                                keepType = 'ALL', VAFratio = 0, SBscore = 0,
                                maxIndelLen = Inf, minInterval = 0, 
                                tagFILTER = NULL, ExAC = FALSE, ESP6500 = FALSE, 
                                gnomAD = FALSE, dbSNP = FALSE)
  # KIRC
  }else if(cancerType == 'KIRC'){
    mafFiltered <- mutFilterCom(maf, tumorDP = 0, normalDP = 0, VAF = 0,
                                dbVAF = 0, dbSNP = TRUE, COSMIConly = FALSE,
                                Genomesprojects1000 = FALSE, tumorAD = 0, 
                                VAFratio = 0, maxIndelLen = Inf, minInterval = 0, 
                                tagFILTER = NULL, ExAC = FALSE, ESP6500 = FALSE, 
                                gnomAD = FALSE, keepType = 'ALL', SBscore = 0,)
  # KIRP
  }else if(cancerType == 'KIRP'){
    mafFiltered <- mutFilterCom(maf, tumorDP = 8, normalDP = 6, VAF = 0.07,
                                dbVAF = 0, dbSNP = TRUE, COSMIConly = TRUE,
                                Genomesprojects1000 = TRUE, tumorAD = 0, 
                                VAFratio = 0, maxIndelLen = Inf, minInterval = 0, 
                                tagFILTER = NULL, ExAC = TRUE, ESP6500 = FALSE, 
                                gnomAD = FALSE, keepType = 'ALL', SBscore = 0)
  # LCML
  }else if(cancerType == 'LCML'){
    mafFiltered <- mutFilterCom(maf, tumorDP = 0, normalDP = 0, VAF = 0.2,
                                dbVAF = 0, dbSNP = FALSE, COSMIConly = FALSE,
                                Genomesprojects1000 = TRUE, tumorAD = 0, 
                                VAFratio = 0, maxIndelLen = Inf, minInterval = 0, 
                                tagFILTER = NULL, ExAC = FALSE, ESP6500 = FALSE, 
                                gnomAD = FALSE, keepType = 'ALL', SBscore = 0)
  # LAML
  }else if(cancerType == 'LAML'){
    mafFiltered <- mutFilterCom(maf, tumorDP = 0, normalDP = 0, VAF = 0,
                                dbVAF = 0, dbSNP = TRUE, COSMIConly = FALSE,
                                Genomesprojects1000 = FALSE, tumorAD = 3, 
                                VAFratio = 0, maxIndelLen = Inf, minInterval = 0, 
                                ExAC = TRUE, ESP6500 = TRUE, 
                                gnomAD = FALSE, keepType = 'ALL', SBscore = 0)
  
  # LIHC
  }else if(cancerType == 'LIHC'){
    mafFiltered <- mutFilterCom(maf, tumorDP = 15, normalDP = 15, VAF = 0.1,
                                dbVAF = 0, dbSNP = TRUE, COSMIConly = TRUE,
                                Genomesprojects1000 = TRUE, tumorAD = 0, 
                                VAFratio = 0, maxIndelLen = Inf, minInterval = 0, 
                                tagFILTER = NULL, ExAC = FALSE, ESP6500 = FALSE, 
                                gnomAD = FALSE, keepType = 'ALL', SBscore = 0)
  }else{
    # turn off the switch for report
    withType <- FALSE
    stop('Invaild cancer type detected, please input a vaild cancer type.')
  }
  return(mafFiltered)
}
