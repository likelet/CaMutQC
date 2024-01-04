#' vcfToMAF
#' @description Format transformation from VCF to MAF.
#'
#' @param vcfFile Directory of a VCF file, or the path to several VCF files
#' that is going to be transformed. Files should be in .vcf or .vcf.gz format.
#' @param multiVCF Logical, whether the input is a path that leads to several
#' VCFs that come from multi-region/sample/caller sequencing. Default: FALSE
#' @param inputStrelka The type of variants ('indel' or 'SNV') in VCF file 
#' if it is from Strelka. Default: FALSE
#' @param writeFile Whether to directly write MAF file to the disk. If FALSE,
#' a MAF data frame will be returned. If TRUE, a MAF file will be saved.
#' Default: FALSE.
#' @param MAFfile File name of the exported MAF file, if writeFile is set as
#' TRUE.
#' @param MAFdir Directory of the exported MAF file, if writeFile is set as
#' TRUE.
#' @param tumorSampleName Name of the tumor sample(s) in the VCF file(s).
#' If it is set as 'Extracted', tumorSampleName would be extracted 
#' automatically from the VCF file. Default: 'Extracted'.
#' @param normalSampleName Name the normal sample in the VCF file.
#' If it is set as 'Extracted', normalSampleName would be extracted 
#' automatically from the VCF file. Default: 'Extracted'.
#' @param ncbiBuild The reference genome used for the alignment, which will be
#' presented as value in 'NCBIbuild' column in MAF file. Default: 'GRCh38'.
#' @param MAFcenter One or more genome sequencing center reporting the variant,
#' which will be presented as value in 'Center' column in MAF. Default: '.'.
#' @param MAFstrand Genomic strand of the reported allele, which will be
#' presented as value in 'Strand' column in MAF file. Default: '+'.
#' @param filterGene Logical. Whether to filter variants without Hugo Symbol.
#' Default: FALSE
#' @param simplified Logical. Whether to extract the first thirteen columns 
#' after converting to MAF file. Default: FALSE
#'
#' @import vcfR org.Hs.eg.db clusterProfiler stringr dplyr utils
#' @importFrom stats na.omit
#' @return A detailed MAF data frame
#' @export vcfToMAF
#'
#' @examples
#' maf <- vcfToMAF(system.file("extdata", "GC48-2_mutect2.vep.vcf",
#' package = "CaMutQC"))


vcfToMAF <- function(vcfFile, multiVCF = FALSE, inputStrelka = FALSE,
                     writeFile = FALSE, MAFfile = 'MAF.maf', MAFdir = './',
                     tumorSampleName = 'Extracted',
                     normalSampleName = 'Extracted', ncbiBuild = 'Extracted',
                     MAFcenter = '.', MAFstrand = '+', filterGene = FALSE,
                     simplified = FALSE){
    # if inputs are multi-sample/multi-caller VCFs
    if (multiVCF){
      filenames <- list.files(vcfFile, pattern="*.vep.vcf", full.names=TRUE)
      if (length(filenames) == 0){
        stop('No VCF file detected under the path you offered.')
      }else{
        # transform each vcf to maf separately
        dfs <- lapply(filenames, vcfhelper, tumorSampleName = tumorSampleName,
                      normalSampleName = normalSampleName,ncbiBuild = ncbiBuild,
                      MAFcenter = MAFcenter, MAFstrand = MAFstrand)
        maf <- do.call("rbind", dfs)
      }
    }else{
      # read vcf File
      ## check whether the file exists
      if (!(file_test('-f', vcfFile))){
        stop('VCF file doesn\'t exist, please check your input.')
      }
      ## check the file format
      nameSplit <- strsplit(vcfFile, split = ".", fixed = TRUE)[[1]]
      suffixLast <- nameSplit[length(nameSplit)]
      suffixSecond <- nameSplit[length(nameSplit) - 1]
      if (suffixLast != 'vcf') {
        stop('Please input file in .vcf or .gz.vcf format')
      }
      maf <- vcfhelper(vcfFile, tumorSampleName = tumorSampleName,
                       normalSampleName=normalSampleName, ncbiBuild = ncbiBuild,
                       MAFcenter = MAFcenter, MAFstrand = MAFstrand,
                       inputStrelka = inputStrelka)
    }
    # simplify MAF file
    if (simplified) {
      maf <- maf[, seq_len(13)]
    }
    if (filterGene){maf <- maf[which(maf$Hugo_Symbol != ''), ]}
    if (writeFile) {
      message('The generated MAF file has been saved.')
      write.table(maf, paste0(MAFdir, MAFfile), sep = "\t",
                  quote = FALSE, row.names = FALSE)
      return(maf)
    }else{return(maf)}
}

## helper function for the transformation from VCF to MAF
vcfhelper <- function(vcfFile, tumorSampleName = 'Extracted',
                      normalSampleName = 'Extracted', ncbiBuild = 'Extracted',
                      MAFcenter = '.', MAFstrand = '+', inputStrelka = FALSE) {
  # message('Loading VCF data...')
  invisible(capture.output(Anno_Vcf <- read.vcfR(vcfFile)))
  mes <- paste0(vcfFile, ' has been loaded successfully!')
  message(mes)

  ## convert to data frames and store useful information
  vcf_header <- as.data.frame(Anno_Vcf@meta)
  vcf_main <- as.data.frame(Anno_Vcf@fix)
  vcf_additional <- as.data.frame(Anno_Vcf@gt)

  # build data frame in MAF format
  maf <- as.data.frame(matrix(ncol = 46, nrow = nrow(vcf_main)))
  colnames(maf) <- c('Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build',
                     'Chromosome', 'Start_Position', 'End_Position', 'Strand',
                     'Variant_Classification', 'Variant_Type',
                     'Reference_Allele','Tumor_Seq_Allele1','Tumor_Seq_Allele2',
                     'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode',
                     'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1',
                     'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1',
                     'Tumor_Validation_Allele2','Match_Norm_Validation_Allele1',
                     'Match_Norm_Validation_Allele2', 'Verification_Status',
                     'Validation_Status', 'Mutation_Status', 'Sequencing_Phase',
                     'Sequence_Source', 'Validation_Method', 'Score',
                     'BAM_File', 'Sequencer', 'Tumor_Sample_UUID',
                     'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp','HGVSp_Short',
                     'Transcript_ID', 'Exon_Number', 't_depth', 't_ref_count',
                     't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count',
                     'all_effects'
  )

  message('VCF to MAF conversion is in process...')
  # assign chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele1 info
  maf[, 5] <- vcf_main$CHROM
  maf[, 6] <- vcf_main$POS
  maf[, 11] <- vcf_main$REF
  maf[, 13] <- vcf_main$ALT

  # get INFO information
  INFOall <- readINFO(vcf_header, vcf_main)
  INFOFrame <- as.data.frame(INFOall[[1]])
  INFOinfo <- as.data.frame(INFOall[[2]])

  # process CSQ column
  ## extract CSQ header and set as colnames
  CSQ_headers <- strsplit(INFOFrame[which(INFOFrame$ID == 'CSQ'), 4],
                          split = "Format: ")[[1]][2]
  CSQ_headers <- strsplit(CSQ_headers, split = "\\|")[[1]]
  CSQ_info <- as.data.frame(matrix(ncol =length(CSQ_headers), nrow = nrow(maf)))
  colnames(CSQ_info) <- CSQ_headers

  ## use built functions to select proper transcripts
  CSQ_info <- selectTrans(maf, CSQ_info, INFOinfo)

  # fill in other information in MAF file
  ## ID conversion
  IDs <- bitr(CSQ_info$Gene, fromType = "ENSEMBL", toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
  ## center
  maf[, 3] <- MAFcenter
  ## NCBI build
  if (ncbiBuild == 'Extracted') {
    maf[, 4] <- na.omit(str_extract(vcf_header$`Anno_Vcf@meta`,
                                    pattern = "(GRCh)[:digit:]+"))[1]
  } else {
    maf[, 4] <- ncbiBuild
  }

  ## tumor sample name and normal sample name
  if (tumorSampleName == 'Extracted'){
    # mutect format: tumor_sample=
    tumorSampleLine <- vcf_header$`Anno_Vcf@meta`[grep('tumor_sample=',
                                                    vcf_header$`Anno_Vcf@meta`)]
    if (length(tumorSampleLine) != 0){
      tumorSampleName <- strsplit(tumorSampleLine, split = '=')[[1]][2]
    }else{
      # MuSE format: TUMOR=
      tumorSampleLine <- vcf_header$`Anno_Vcf@meta`[grep('TUMOR=',
                                                  vcf_header$`Anno_Vcf@meta`)]
      if (length(tumorSampleLine) != 0){
        tumorSampleName <- strsplit(strsplit(tumorSampleLine, 
                                    split = ',')[[1]][1], split = "=")[[1]][3]
      }else{
        tumorSampleName <- 'TUMOR'
      }
    }
  }

  if (normalSampleName == 'Extracted'){
    # mutect format: normal_sample=
    normalSampleLine <- vcf_header$`Anno_Vcf@meta`[grep('normal_sample=',
                                                vcf_header$`Anno_Vcf@meta`)]
    if (length(normalSampleLine) != 0){
      normalSampleName <- strsplit(normalSampleLine, split = '=')[[1]][2]
    }else{
      # MuSE format: NORMAL=
      normalSampleLine <- vcf_header$`Anno_Vcf@meta`[grep('NORMAL=',
                                                  vcf_header$`Anno_Vcf@meta`)]
      if (length(normalSampleLine) != 0){
        normalSampleName <- strsplit(strsplit(normalSampleLine, 
                                      split = ',')[[1]][1], split = "=")[[1]][3]
      }else{
        normalSampleName <- 'NORMAL'
      }
    }
  }
  colnames(vcf_additional)[which(colnames(vcf_additional) == 'TUMOR')] <- tumorSampleName
  colnames(vcf_additional)[which(colnames(vcf_additional) == 'NORMAL')] <- normalSampleName
  ## strand
  maf[, 8] <- MAFstrand
  ## Hugo_Symbol
  maf[, 1] <- CSQ_info[, 4]
  ## add VAF column
  maf <- cbind(maf, VAF = 0)
  maf[, c(15, 18:34, 37, 46)] <- '.'
  for (i in seq_len(nrow(maf))) {
    ## ENTREZID
    maf[i, 2] <- IDs$ENTREZID[which(IDs$ENSEMBL == CSQ_info$Gene[i])][1]
    ## get correct variant Position, Variant_Type, Ref allele and Alt allele
    pos_type <- getVarFeature(maf[i, 6], maf[i, 11], maf[i, 13], CSQ_info[i, 1])
    maf[i, 6] <- pos_type[[1]] # start
    maf[i, 7] <- pos_type[[2]] # end
    maf[i, 10] <- pos_type[[3]] # Variant_Type
    maf[i, 11] <- pos_type[[4]] # correct ref allele
    maf[i, 13] <- pos_type[[5]] # correct alt allele
    inframe <- pos_type[[6]] # inframe indicator
    ## Variant_Class
    ### select consequence in CSQ first
    cons <- strsplit(CSQ_info[i, 2], split = '&')[[1]]
    CSQ_info[i, 2] <- cons[which.min(unlist(vapply(cons,
                                                  GetConsequencePriority, 10)))]
    maf[i, 9] <- getVarClass(CSQ_info[i, 2], maf[i, 10], inframe)
    ## get AD and DP in INFO or FORMAT
    ## if RD field is contained in the INFO column
    if (any(strsplit(vcf_additional[i, 1], ":")[[1]] == 'RD')) {
      RD_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'RD'
      AD_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'AD'
      DP_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'DP'
      maf[i, "t_ref_count"] <- strsplit(strsplit(vcf_additional[i, tumorSampleName],
                                                ":")[[1]][RD_loc], ",")[[1]]
      maf[i, "t_alt_count"] <- strsplit(strsplit(vcf_additional[i, tumorSampleName],
                                                 ":")[[1]][AD_loc], ",")[[1]]
      maf[i, "n_ref_count"] <- strsplit(strsplit(vcf_additional[i, normalSampleName],
                                                 ":")[[1]][RD_loc], ",")[[1]]
      maf[i, "n_alt_count"] <- strsplit(strsplit(vcf_additional[i, normalSampleName],
                                                 ":")[[1]][AD_loc], ",")[[1]]
      maf[i, "t_depth"] <- strsplit(strsplit(vcf_additional[i, tumorSampleName],
                                             ":")[[1]][DP_loc], ",")[[1]]
      maf[i, "n_depth"] <- strsplit(strsplit(vcf_additional[i, normalSampleName],
                                            ":")[[1]][DP_loc], ",")[[1]]

      maf[i, 'VAF'] <- as.numeric(maf[i, "t_alt_count"])/
        as.numeric(maf[i, "t_depth"])
    }else{
      AD_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'AD'
      ## if AD field contains information from both ref and alt
      if (any(AD_loc)){
        AD <- strsplit(strsplit(vcf_additional[i, tumorSampleName],
                                ":")[[1]][AD_loc], ",")[[1]][2]
        if (length(grep('DP', strsplit(vcf_additional[i, 1], ":")[[1]]))){
          DP_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'DP'
          DP <- as.numeric(strsplit(vcf_additional[i, tumorSampleName],
                                    ":")[[1]][DP_loc])
          tDP <- as.numeric(strsplit(vcf_additional[i, tumorSampleName],
                                     ":")[[1]][DP_loc])
          nDP <- as.numeric(strsplit(vcf_additional[i, normalSampleName],
                                     ":")[[1]][DP_loc])
        }else if('DP' %in% colnames(INFOinfo)){
          DP <- as.numeric(INFOinfo[i, 'DP'])
          tDP <- DP
          nDP <- sum(as.numeric(strsplit(strsplit(vcf_additional[i, normalSampleName],
                                                  ":")[[1]][AD_loc], ",")[[1]]))
        }
        ## t_depth, n_depth, t_ref_count, t_alt_count, n_ref_count, n_alt_count
        tRefAD <- as.numeric(strsplit(strsplit(vcf_additional[i, tumorSampleName],
                                               ":")[[1]][AD_loc], ",")[[1]][1])
        nRefAD <- as.numeric(strsplit(strsplit(vcf_additional[i, normalSampleName],
                                               ":")[[1]][AD_loc], ",")[[1]][1])
        nAltAD <- as.numeric(strsplit(strsplit(vcf_additional[i, normalSampleName],
                                               ":")[[1]][AD_loc], ",")[[1]][2])
        maf[i, "t_ref_count"] <- tRefAD
        maf[i, "t_alt_count"] <- as.numeric(AD)
        maf[i, "n_ref_count"] <- nRefAD
        maf[i, "n_alt_count"] <- nAltAD
        if (!(exists('tDP'))){
          maf[i, "t_depth"] <- tRefAD + as.numeric(AD)
          maf[i, "n_depth"] <- nRefAD + nAltAD
          maf[i, 'VAF'] <- as.numeric(AD)/maf[i, "t_depth"]
        }else{
          maf[i, "t_depth"] <- tDP
          maf[i, "n_depth"] <- nDP
          maf[i, 'VAF'] <- as.numeric(AD)/DP
        }
      }else{
        # if no AD field contained in the VCF, it should be output from Strelka
        if (inputStrelka == 'indel'){
          TIR_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'TIR'
          TAR_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'TAR'
          tTIR <- strsplit(vcf_additional[i, tumorSampleName],":")[[1]][TIR_loc]
          tTAR <- strsplit(vcf_additional[i, tumorSampleName],":")[[1]][TAR_loc]
          nTIR <- strsplit(vcf_additional[i, normalSampleName],":")[[1]][TIR_loc]
          nTAR <- strsplit(vcf_additional[i, normalSampleName],":")[[1]][TAR_loc]
          # use tire1 information
          tAltAD <- as.numeric(strsplit(tTIR, ",")[[1]][1])
          tRefAD <- as.numeric(strsplit(tTAR, ",")[[1]][1])
          nAltAD <- as.numeric(strsplit(nTIR, ",")[[1]][1])
          nRefAD <- as.numeric(strsplit(nTAR, ",")[[1]][1])
        }else if (inputStrelka == 'SNV'){
          refAllele <- paste0(vcf_main[i, 'REF'], "U")
          altAllele <- paste0(vcf_main[i, 'ALT'], "U")
          ref_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == refAllele
          alt_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == altAllele
          tRef <- strsplit(vcf_additional[i, tumorSampleName],":")[[1]][ref_loc]
          tAlt <- strsplit(vcf_additional[i, tumorSampleName],":")[[1]][alt_loc]
          nRef <- strsplit(vcf_additional[i, normalSampleName],":")[[1]][ref_loc]
          nAlt <- strsplit(vcf_additional[i, normalSampleName],":")[[1]][alt_loc]
          tAltAD <- as.numeric(strsplit(tAlt, ",")[[1]][1])
          tRefAD <- as.numeric(strsplit(tRef, ",")[[1]][1])
          nAltAD <- as.numeric(strsplit(nAlt, ",")[[1]][1])
          nRefAD <- as.numeric(strsplit(nRef, ",")[[1]][1])
        }else{
          mes <- paste0("If your input comes from Strelka,",
                        " set 'indel' or 'SNV' for 'inputStrelka' parameter")
          stop(mes)
        }
        maf[i, "t_ref_count"] <- tRefAD
        maf[i, "t_alt_count"] <- tAltAD
        maf[i, "n_ref_count"] <- nRefAD
        maf[i, "n_alt_count"] <- nAltAD
        maf[i, "t_depth"] <- tRefAD + tAltAD
        maf[i, "n_depth"] <- nRefAD + nAltAD
        maf[i, 'VAF'] <- tAltAD/maf[i, "t_depth"]
      }
    }
    ## fill in dbSNP_RS
    if(nchar(CSQ_info[i, 'Existing_variation']) == 0) {
      maf[i, 'dbSNP_RS'] <- 'novel'
    } else if (str_detect(CSQ_info[i, 'Existing_variation'], 'rs')) {
      extVar <- strsplit(CSQ_info[i, 'Existing_variation'], "&")[[1]]
      dbSNP_loc <- str_detect(extVar, 'rs')
      maf[i, 'dbSNP_RS'] <- paste(extVar[dbSNP_loc], collapse = ",")
    } else {
      maf[i, 'dbSNP_RS'] <- '.'
    }
    ## HGVSc
    if (nchar(CSQ_info[i, 'HGVSc']) == 0) {
      maf[i, 'HGVSc'] <- ''
    }else{
      maf[i, 'HGVSc'] <- strsplit(CSQ_info[i, 'HGVSc'], split = ":")[[1]][2]
    }
    ## HGVSp
    if (nchar(CSQ_info[i, 'HGVSp']) == 0) {
      maf[i, 'HGVSp'] <- ''
    }else{
      maf[i, 'HGVSp'] <- strsplit(CSQ_info[i, 'HGVSp'], split = ":")[[1]][2]
    }
  }
  ## Transcript_ID
  maf[, 'Transcript_ID'] <- CSQ_info[, 'Feature']
  ## Exon_Number
  maf[, 'Exon_Number'] <- CSQ_info[, 'EXON']
  ## set Tumor_Seq_Allele1 same as ref
  maf[, 12] <- maf[, 11]
  maf[, 'Tumor_Sample_Barcode'] <- tumorSampleName
  maf[, 'Matched_Norm_Sample_Barcode'] <- normalSampleName
  maf$t_alt_count <- as.numeric(maf$t_alt_count)
  maf$t_ref_count <- as.numeric(maf$t_ref_count)
  maf$t_depth <- as.numeric(maf$t_depth)
  maf$n_alt_count <- as.numeric(maf$n_alt_count)
  maf$n_ref_count <- as.numeric(maf$n_ref_count)
  maf$n_depth <- as.numeric(maf$n_depth)
  maf1 <- jointMAF(maf[ ,seq_len(46)], CSQ_info, vcf_main)
  # remane NORMAL and TUMOR column so that vcf from different callers can merge
  colnames(vcf_additional)[which(colnames(vcf_additional) 
                                 == tumorSampleName)] <- 'tumorSampleInfo'
  colnames(vcf_additional)[which(colnames(vcf_additional) 
                                 == normalSampleName)] <- 'normalSampleInfo'
  maf <- cbind(maf1, VAF = maf[ ,'VAF'], vcf_additional)
  # change column type
  maf$Start_Position <- as.numeric(maf$Start_Position)
  maf$End_Position <- as.numeric(maf$End_Position)
  # change column type without displaying warnings
  maf <- maf %>% mutate_at(vars(78:86, 101:117), as.numeric)
  maf <- cbind(maf, CaTag = '0')
  rownames(maf) <- seq_len(nrow(maf))
  message('VCF to MAF conversion has been done successfully!')
  return(maf)
}
