#' vcfToMAF
#' @description Format transformation from VCF to MAF.
#'
#' @param vcfFile Directory of a VCF file, or the path to several VCF files
#' that is going to be transformed. Files should be in .vcf or .vcf.gz format.
#' @param multiSample Logical, whether the input is a path that leads to several
#' VCFs that come from multi-region/sample sequencing. Default: FALSE
#' @param writeFile Whether to directly write MAF file to the disk. If FALSE,
#' a MAF data frame will be returned. If TRUE,
#' a MAF file will be saved in your disk. Default: FALSE.
#' @param MAFfile File name of the exported MAF file, if writeFile is set as
#' TRUE.
#' @param MAFdir Directory of the exported MAF file, if writeFile is set as
#' TRUE.
#' @param tumorSampleName Name of the tumor sample(s) in the VCF file(s).
#' If it is set as 'Extracted', tumorSampleName would be extracted automatically
#' from the VCF file. Default: 'Extracted'.
#' @param normalSampleName Name the normal sample in the VCF file.
#' If it is set as 'Extracted', normalSampleName would be extracted automatically
#' from the VCF file. Default: 'Extracted'.
#' @param ncbiBuild The reference genome used for the alignment, which will be
#' presented as value in 'NCBIbuild' column in MAF file. Default: 'GRCh38'.
#' @param MAFcenter One or more genome sequencing center reporting the variant,
#' which will be presented as value in 'Center' column in MAF file. Default: '.'.
#' @param MAFstrand Genomic strand of the reported allele, which will be
#' presented as value in 'Strand' column in MAF file. Default: '+'.
#' @param filterGene Logical. Whether to filter variants without Hugo Symbol.
#' Default: TRUE
#' @param simplified Logical. Whether to extract the first thirteen columns after
#' converting to MAF file. Default: FALSE
#'
#' @import vcfR org.Hs.eg.db clusterProfiler stringr dplyr
#' @return A MAF data frame
#' @export vcfToMAF
#'
#' @examples
#' maf <- vcfToMAF(system.file("extdata", "GC48-2_mutect2.vep.vcf",
#' package = "CaMutQC"))


vcfToMAF <- function(vcfFile, multiSample = FALSE, writeFile = FALSE,
                     MAFfile = 'MAF.maf', MAFdir = './',
                     tumorSampleName = 'Extracted',
                     normalSampleName = 'Extracted', ncbiBuild = 'Extracted',
                     MAFcenter = '.', MAFstrand = '+', filterGene = TRUE,
                     simplified = FALSE){

  # if inputs are multi-sample VCFs
  if (multiSample){
    filenames <- list.files(vcfFile, pattern="*.vep.vcf", full.names=TRUE)
    if (length(filenames) == 0){
      stop('No VCF file detected under the path you offered.')
    }else{
      dfs <- lapply(filenames, vcfhelper, tumorSampleName = tumorSampleName,
                    normalSampleName = normalSampleName, ncbiBuild = ncbiBuild,
                    MAFcenter = MAFcenter, MAFstrand = MAFstrand)
      #if (length(dfs) == length(tumorSampleNames)){
       # names(dfs) <- tumorSampleNames
      #}else{
       # warning('Length of objects doesn\'t match')
      #}
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
                     normalSampleName = normalSampleName, ncbiBuild = ncbiBuild,
                     MAFcenter = MAFcenter, MAFstrand = MAFstrand)
  }

  # simplify MAF file
  if (simplified) {
    maf <- maf[, 1:13]
    warning(paste0('Simplifying MAF file will only keep first thirteen columns, ',
                   'so following filtration cannot be proceeded on simplified MAF file.'))
  }

  if (filterGene){
    maf <- maf[which(maf$Hugo_Symbol != ''), ]
  }

  if (writeFile) {
    message('The generated MAF file has been saved.')
    write.table(maf, paste0(MAFdir, MAFfile), sep = "\t",
                quote = FALSE, row.names = FALSE)

  } else{
    return(maf)
  }

}



vcfhelper <- function(vcfFile, tumorSampleName = 'Extracted',
                      normalSampleName = 'Extracted', ncbiBuild = 'Extracted',
                      MAFcenter = '.', MAFstrand = '+') {
  message('Loading VCF data...')
  invisible(capture.output(Anno_Vcf <- read.vcfR(vcfFile)))
  message('VCF data has been loaded successfully!')

  ## convert to data frames and store useful information
  vcf_header <- as.data.frame(Anno_Vcf@meta)
  vcf_main <- as.data.frame(Anno_Vcf@fix)
  vcf_additional <- as.data.frame(Anno_Vcf@gt)

  # build data frame in MAF format
  maf <- as.data.frame(matrix(ncol = 46, nrow = nrow(vcf_main)))
  colnames(maf) <- c('Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build',
                     'Chromosome', 'Start_Position', 'End_Position', 'Strand',
                     'Variant_Classification', 'Variant_Type',
                     'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2',
                     'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode',
                     'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1',
                     'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1',
                     'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1',
                     'Match_Norm_Validation_Allele2', 'Verification_Status',
                     'Validation_Status', 'Mutation_Status', 'Sequencing_Phase',
                     'Sequence_Source', 'Validation_Method', 'Score',
                     'BAM_File', 'Sequencer', 'Tumor_Sample_UUID',
                     'Matched_Norm_Sample_UUID', 'HGVSc', 'HGVSp', 'HGVSp_Short',
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
  INFOFrame <- as.data.frame(readINFO(vcf_header, vcf_main)[[1]])
  INFOinfo <- as.data.frame(readINFO(vcf_header, vcf_main)[[2]])

  # process CSQ column
  ## extract CSQ header and set as colnames
  CSQ_headers <- strsplit(INFOFrame[which(INFOFrame$ID == 'CSQ'), 4],
                          split = "Format: ")[[1]][2]
  CSQ_headers <- strsplit(CSQ_headers, split = "\\|")[[1]]
  CSQ_info <- as.data.frame(matrix(ncol = length(CSQ_headers), nrow = nrow(maf)))
  colnames(CSQ_info) <- CSQ_headers

  ## use built functions to select proper transcripts
  CSQ_info <- selectTrans(maf, CSQ_info, INFOinfo)

  # fill in other information in MAF file
  ## ID conversion
  IDs <- suppressWarnings(suppressMessages(bitr(CSQ_info$Gene,
                                                fromType = "ENSEMBL",
                                                toType = "ENTREZID",
                                                OrgDb = org.Hs.eg.db)))

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
    tumorSampleLine <- vcf_header$`Anno_Vcf@meta`[grep('tumor_sample=',
                                                       vcf_header$`Anno_Vcf@meta`)]

    if (length(tumorSampleLine) == 0){
      #('Please check the input VCF file to make sure that it contains
           #\'tumor_sample=\' in the header.')
      # if no ID was given in the header of VCF file, CaMutQC will use colnamaes
      # of the VCF main section
      tumorSampleName <- 'TUMOR'

    }else{
      tumorSampleName <- strsplit(tumorSampleLine, split = '=')[[1]][2]

    }
  }

  if (normalSampleName == 'Extracted'){
    normalSampleLine <- vcf_header$`Anno_Vcf@meta`[grep('normal_sample=',
                                                        vcf_header$`Anno_Vcf@meta`)]

    if (length(normalSampleLine) == 0){
      #('Please check the input VCF file to make sure that it contains
      #\'normal_sample=\' in the header.')
      # if no ID was given in the header of VCF file, CaMutQC will use colnamaes
      # of the VCF main section
      normalSampleName <- 'NORMAL'
    }else{
      normalSampleName <- strsplit(normalSampleLine, split = '=')[[1]][2]

    }
  }

  ### check the consistence of tumorSampleName/normalSampleName
  ### with colnames of FORMAT column
  if (!(tumorSampleName %in% colnames(vcf_additional))){
    stop('Tumor Sample Name is invalid.
     Set as \'Extracted\' if you want tumorSampleName
         to be extracted from VCF file')
  }else if(!(normalSampleName %in% colnames(vcf_additional))){
    stop('Normal Sample Name is invalid.
     Set as \'Extracted\' if you want normalSampleName
         to be extracted from VCF file')
  }

  ## strand
  maf[, 8] <- MAFstrand

  ## Hugo_Symbol
  maf[, 1] <- CSQ_info[, 4]

  ## add VAF column
  maf <- cbind(maf, VAF = 0)

  ## Tumor_Seq_Allele2
  #if (length(strsplit(maf2[i, 13], split = ",")[[1]]) != 1) {
  #maf2[i, 13] <- selectAlt(maf2[i, ], vcf_additional[i, ],
  #tumorSampleName, normalSampleName)
  #}

  maf[, c(15, 18:34, 37, 46)] <- '.'

  for (i in 1:nrow(maf)) {

    ## ENTREZID
    maf[i, 2] <- IDs$ENTREZID[which(IDs$ENSEMBL == CSQ_info$Gene[i])][1]


    ## endPos and Variant_Type
    maf[i, 7] <- getVariantType(maf[i, 6], maf[i, 11], maf[i, 13])[1]
    maf[i, 10] <- getVariantType(maf[i, 6], maf[i, 11], maf[i, 13])[2]
    inframe <- switch(getVariantType(maf[i, 6], maf[i, 11], maf[i, 13])[3],
                      'TRUE' = 1,
                      'FALSE' = 0)

    ## Variant_Class
    ### select consequence in CSQ first
    cons <- strsplit(CSQ_info[i, 2], split = '&')[[1]]
    CSQ_info[i, 2] <- cons[which.max(unlist(sapply(cons,
                                                  GetConsequencePriority)))]
    #if (length(strsplit(CSQ_info[i, 2], split = '&')[[1]]) > 1) {
      #cons1 <- strsplit(CSQ_info[i, 2], split = '&')[[1]][1]
      #cons2 <- strsplit(CSQ_info[i, 2], split = '&')[[1]][2]
      #if (GetConsequencePriority(cons1) >= GetConsequencePriority(cons2)) {
        #CSQ_info[i, 2] <- cons1
      #} else {
        #CSQ_info[i, 2] <- cons2
      #}
    #}

    maf[i, 9] <- getVariantClassification(CSQ_info[i, 2], maf[i, 10], inframe)

    ## fill in VAF
    #if (length(grep('AF', strsplit(vcf_additional[i, 1], ":")[[1]]))){
    #AF_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'AF'
    #VAFs <- strsplit(vcf_additional[i, tumorSampleName], ":")[[1]][AF_loc]
    #if (length(grep(',', VAFs))){
    #maf[i, 'VAF'] <- as.numeric(strsplit(VAFs, ',')[[1]][1])
    #}else{
    #maf[i, 'VAF'] <- as.numeric(VAFs)
    #}
    #}

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
      ## if AD field contains information from both ref and alt
      AD_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'AD'
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
      #tAltAD <- as.numeric(strsplit(strsplit(vcf_additional[i, tumorSampleName],
                                             #":")[[1]][AD_loc], ",")[[1]][2])
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

  ## Tumor_Seq_Allele1, set Tumor_Seq_Allele1 same as ref
  maf[, 12] <- maf[, 11]

  maf[, 'Tumor_Sample_Barcode'] <- tumorSampleName
  maf[, 'Matched_Norm_Sample_Barcode'] <- normalSampleName


  maf1 <- jointMAF(maf[ ,1:46], CSQ_info, vcf_main)
  maf <- cbind(maf1, VAF = maf[ ,'VAF'], vcf_additional)


  # maf <- cbind(maf, isIndelAround = 0)
  maf$Start_Position <- as.numeric(maf$Start_Position)
  maf$End_Position <- as.numeric(maf$End_Position)
  maf$t_alt_count <- as.numeric(maf$t_alt_count)
  maf$t_ref_count <- as.numeric(maf$t_ref_count)
  maf$t_depth <- as.numeric(maf$t_depth)
  maf$n_alt_count <- as.numeric(maf$n_alt_count)
  maf$n_ref_count <- as.numeric(maf$n_ref_count)
  maf$n_depth <- as.numeric(maf$n_depth)

  colnames(maf)[which(colnames(maf) == tumorSampleName)] <- 'tumorSampleInfo'
  maf <- cbind(maf, CaTag = '0')
  rownames(maf) <- 1:nrow(maf)
  message('VCF to MAF conversion has been done successfully!')
  return(maf)
}
