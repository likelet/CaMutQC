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
#' @param dbSNP Whether to filter variants listed in dbSNP. Default: FALSE.
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
#' @examples
#' maf <- vcfToMAF(system.file("extdata/Multi-sample",
#' "SRR3670028.somatic.filter.HC.vep.vcf",package = "CaMutQC"))
#' mafF <- mutFilterCan(maf, cancerType = 'BRCA')

mutFilterCan <- function(maf, cancerType, tumorDP = 0, normalDP = 0,
                         tumorAD = 0, VAF = 0, VAFratio = 0,
                         tumorSampleName = 'Extracted', SBmethod = 'SOR',
                         SBscore = Inf, maxIndelLen = Inf, minInterval = 0,
                         tagFILTER = NULL, dbVAF = 0, ExAC = FALSE,
                         Genomesprojects1000 = FALSE, ESP6500 = FALSE,
                         gnomAD = FALSE, dbSNP = FALSE, COSMIConly = FALSE,
                         keepType = 'ALL', bedFile = NULL, bedFilter = TRUE,
                         mutFilter = FALSE, selectCols = FALSE, report = TRUE,
                         reportFile = 'FilterReport.html', reportDir = './',
                         TMB = FALSE) {

  # turn on the switch for report
  # withType <- TRUE
  # BLCA
  if (cancerType == 'BLCA'){
    mafFiltered <- mutFilterCom(maf, SBmethod = 'Fisher', SBscore = 20,
                                minInterval = 30, tumorDP = 10,
                                Genomesprojects1000 = TRUE, ExAC = TRUE,
                                normalDP = 10, tumorAD = 5, VAF = VAF,
                                VAFratio = VAFratio, maxIndelLen = maxIndelLen,
                                tumorSampleName = tumorSampleName,
                                tagFILTER = tagFILTER, dbVAF = dbVAF,
                                ESP6500 = ESP6500, gnomAD = gnomAD,
                                dbSNP = dbSNP, COSMIConly = COSMIConly,
                                keepType = keepType, bedFile = bedFile,
                                bedFilter = bedFilter, mutFilter = mutFilter,
                                selectCols = selectCols, report = report,
                                reportFile = reportFile, reportDir = reportDir,
                                TMB = TMB)
  # BRCA
  }else if(cancerType == 'BRCA'){
    mafFiltered <- mutFilterCom(maf, tumorAD = 5, VAF = 0.1,
                                dbSNP = TRUE, Genomesprojects1000 = TRUE,
                                ESP6500 = TRUE, COSMIConly = TRUE,
                                tumorDP = 6, normalDP = 6, VAFratio = VAFratio,
                                tumorSampleName = tumorSampleName,
                                SBmethod = SBmethod, ExAC = ExAC,
                                SBscore = SBscore, maxIndelLen = maxIndelLen,
                                minInterval = minInterval, bedFile = bedFile,
                                tagFILTER = tagFILTER, dbVAF = dbVAF,
                                gnomAD = gnomAD, keepType = keepType,
                                bedFilter = bedFilter, mutFilter = mutFilter,
                                selectCols = selectCols, report = report,
                                reportFile = reportFile, reportDir = reportDir,
                                TMB = TMB)
  # COADREAD
  }else if(cancerType == 'COADREAD'){
    mafFiltered <- mutFilterCom(maf, tumorDP = 5, VAF = 0.2,
                                dbSNP = TRUE, Genomesprojects1000 = TRUE,
                                normalDP = normalDP, tumorAD = tumorAD,
                                VAFratio = VAFratio, SBmethod = SBmethod,
                                tumorSampleName = tumorSampleName,
                                SBscore = SBscore, maxIndelLen = maxIndelLen,
                                minInterval = minInterval, ExAC = ExAC,
                                tagFILTER = tagFILTER, dbVAF = dbVAF,
                                ESP6500 = ESP6500, gnomAD = gnomAD,
                                COSMIConly = COSMIConly, keepType = keepType,
                                bedFile = bedFile, bedFilter = bedFilter,
                                mutFilter = mutFilter, selectCols = selectCols,
                                report = report, reportFile = reportFile,
                                reportDir = reportDir, TMB = TMB)
  # UCEC
  }else if(cancerType == 'UCEC'){
    mafFiltered <- mutFilterCom(maf, dbSNP = TRUE, Genomesprojects1000 = TRUE)

  # UCS
  }else if(cancerType == 'UCS'){
    mafFiltered <- mutFilterCom(maf, tumorAD = 5, tumorDP = 12, normalDP = 5,
                                COSMIConly = TRUE)
  # KIRC
  }else if(cancerType == 'KIRC'){
    mafFiltered <- mutFilterCom(maf, dbSNP = TRUE)
  # KIRP
  }else if(cancerType == 'KIRP'){
    mafFiltered <- mutFilterCom(maf, tumorDP = 8, normalDP = 6, VAF = 0.07,
                                dbSNP = TRUE, COSMIConly = TRUE,
                                Genomesprojects1000 = TRUE, ExAC = TRUE)
  # LCML
  }else if(cancerType == 'LCML'){
    mafFiltered <- mutFilterCom(maf, VAF = 0.2, Genomesprojects1000 = TRUE)
  # LAML
  }else if(cancerType == 'LAML'){
    mafFiltered <- mutFilterCom(maf, dbSNP = TRUE, tumorAD = 3,
                                ExAC = TRUE, ESP6500 = TRUE, tagFILTER = 'PASS')

  # LIHC
  }else if(cancerType == 'LIHC'){
    mafFiltered <- mutFilterCom(maf, tumorDP = 15, normalDP = 15, VAF = 0.1,
                                dbSNP = TRUE, COSMIConly = TRUE,
                                Genomesprojects1000 = TRUE)
  }else{
    # turn off the switch for report
    # withType <- FALSE
    stop('Invaild cancer type detected, please input a vaild cancer type.')
  }
  return(mafFiltered)
}
