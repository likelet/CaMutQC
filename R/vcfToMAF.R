#' vcfToMAF
#' @description Format transformation from VCF to MAF.
#'
#' @param vcfFile Directory of the VCF file that is going to be transformed. 
#' It should be in .vcf or .vcf.gz format.
#' @param writeFile Whether to directly write MAF file to the disk. If FALSE,
#' a MAF data frame will be returned. If TRUE, 
#' a MAF file will be saved in your disk. Default: FALSE.
#' @param MAFfile File name of the exported MAF file, if writeFile is set as 
#' TRUE.
#' @param MAFdir Directory of the exported MAF file, if writeFile is set as 
#' TRUE.
#' @param tumorSampleName Name of the tumor sample in the VCF file. 
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
#' 
#' @import vcfR, org.Hs.eg.db, clusterProfiler, stringr, dplyr
#' @return A MAF data frame


vcfToMAF <- function(vcfFile, writeFile = FALSE, MAFfile = 'MAF.maf', 
                     MAFdir = './', tumorSampleName = 'Extracted', 
                     normalSampleName = 'Extracted', ncbiBuild = 'GRCh38', 
                     MAFcenter = '.', MAFstrand = '+'){
  
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
  
  Anno_Vcf <- read.vcfR(vcfFile)
  
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
  IDs <- bitr(CSQ_info$Gene, fromType = "ENSEMBL", 
              toType = "ENTREZID", 
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
    tumorSampleLine <- vcf_header$`Anno_Vcf@meta`[grep('tumor_sample=', 
                                                       vcf_header$`Anno_Vcf@meta`)]
    
    if (length(tumorSampleLine) == 0){
      stop('Please check the input VCF file to make sure that it contains
           \'tumor_sample=\' in the header.')
    }else{
      tumorSampleName <- strsplit(tumorSampleLine, split = '=')[[1]][2]
      
    }
  }
  
  if (normalSampleName == 'Extracted'){
    normalSampleLine <- vcf_header$`Anno_Vcf@meta`[grep('normal_sample=', 
                                                        vcf_header$`Anno_Vcf@meta`)]
    
    if (length(normalSampleLine) == 0){
      stop('Please check the input VCF file to make sure that it contains
           \'normal_sample=\' in the header.')
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
    if (length(strsplit(CSQ_info[i, 2], split = '&')[[1]]) > 1) {
      cons1 <- strsplit(CSQ_info[i, 2], split = '&')[[1]][1]
      cons2 <- strsplit(CSQ_info[i, 2], split = '&')[[1]][2]
      if (GetConsequencePriority(cons1) >= GetConsequencePriority(cons2)) {
        CSQ_info[i, 2] <- cons1
      } else {
        CSQ_info[i, 2] <- cons2
      }
    }
    
    maf[i, 9] <- getVariantClassification(CSQ_info[i, 2], maf[i, 10], inframe)
    
    ## fill in VAF
    DP_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'DP'
    AD_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'AD'
    AD <- strsplit(strsplit(vcf_additional[i, tumorSampleName], ":")[[1]][AD_loc], 
                   ",")[[1]][2]
    DP <- strsplit(vcf_additional[i, tumorSampleName], ":")[[1]][DP_loc]
    maf[i, 'VAF'] <- as.numeric(AD)/as.numeric(DP)
    
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
    
    ## HGVSp_short
    
    
    ## Transcript_ID
    maf[i, 'Transcript_ID'] <- CSQ_info[i, 'Feature']
    
    ## Exon_Number
    
    maf[i, 'Exon_Number'] <- CSQ_info[i, 'EXON']
    
    ## t_depth, n_depth, t_ref_count, t_alt_count, n_ref_count, n_alt_count
    DP_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'DP'
    AD_loc <- strsplit(vcf_additional[i, 1], ":")[[1]] == 'AD'
    tRefAD <- as.numeric(strsplit(strsplit(vcf_additional[i, tumorSampleName], 
                                           ":")[[1]][AD_loc], ",")[[1]][1])
    
    tAltAD <- as.numeric(strsplit(strsplit(vcf_additional[i, tumorSampleName], 
                                           ":")[[1]][AD_loc], ",")[[1]][2])
    tDP <- as.numeric(strsplit(vcf_additional[i, tumorSampleName], 
                               ":")[[1]][DP_loc])
    nDP <- as.numeric(strsplit(vcf_additional[i, normalSampleName], 
                               ":")[[1]][DP_loc])
    nRefAD <- as.numeric(strsplit(strsplit(vcf_additional[i, normalSampleName], 
                                           ":")[[1]][AD_loc], ",")[[1]][1])
    nAltAD <- as.numeric(strsplit(strsplit(vcf_additional[i, normalSampleName], 
                                           ":")[[1]][AD_loc], ",")[[1]][2])
    maf[i, "t_depth"] <- tDP
    maf[i, "n_depth"] <- nDP
    maf[i, "t_ref_count"] <- tRefAD
    maf[i, "t_alt_count"] <- tAltAD
    maf[i, "n_ref_count"] <- nRefAD
    maf[i, "n_alt_count"] <- nAltAD
    
    
    ## Tumor_Seq_Allele1, set Tumor_Seq_Allele1 same as ref
    maf[i, 12] <- maf[i, 11]
  }
  
  maf[, 'Tumor_Sample_Barcode'] <- tumorSampleName
  maf[, 'Matched_Norm_Sample_Barcode'] <- normalSampleName
  
  
  maf <- cbind(maf[, 1:46], Allele = CSQ_info$Allele, Gene = CSQ_info$Gene,
               Feature = CSQ_info$Feature, Feature_type = CSQ_info$Feature_type, 
               One_Consequence = CSQ_info$Consequence, 
               Consequence = CSQ_info$Consequence, 
               cDNA_position = CSQ_info$cDNA_position,
               CDS_position = CSQ_info$CDS_position, 
               Protein_position = CSQ_info$Protein_position,
               Amino_acids = CSQ_info$Amino_acids, Codons = CSQ_info$Codons,
               Existing_variation = CSQ_info$Existing_variation,
               ALLELE_NUM = CSQ_info$ALLELE_NUM, DISTANCE = CSQ_info$DISTANCE,
               TRANSCRIPT_STRAND = CSQ_info$STRAND, SYMBOL = CSQ_info$SYMBOL,
               SYMBOL_SOURCE = CSQ_info$SYMBOL_SOURCE, 
               HGNC_ID = CSQ_info$HGNC_ID, BIOTYPE = CSQ_info$BIOTYPE, 
               CANONICAL = CSQ_info$CANONICAL, CCDS = CSQ_info$CCDS, 
               ENSP = CSQ_info$ENSP, SWISSPROT = CSQ_info$SWISSPROT,
               TREMBL = CSQ_info$TREMBL, UNIPARC = CSQ_info$UNIPARC,
               RefSeq = CSQ_info$RefSeq, SIFT = CSQ_info$SIFT, 
               PolyPhen = CSQ_info$PolyPhen, Exon = CSQ_info$EXON, 
               INTRON = CSQ_info$INTRON, DOMAINS = CSQ_info$DOMAINS,
               GMAF = CSQ_info$gnomAD_AF, AFR_MAF = CSQ_info$AFR_AF, 
               AMR_MAF = CSQ_info$AMR_AF, ASN_MAF = '.', 
               EAS_MAF = CSQ_info$EAS_AF, EUR_MAF = CSQ_info$EUR_AF,
               SAS_MAF = CSQ_info$SAS_AF, AA_MAF = CSQ_info$AA_AF, 
               EA_MAF = CSQ_info$EA_AF, CLIN_SIG = CSQ_info$CLIN_SIG, 
               SOMATIC = CSQ_info$SOMATIC, PUBMED = CSQ_info$PUBMED,
               MOTIF_NAME = CSQ_info$MOTIF_NAME, MOTIF_POS = CSQ_info$MOTIF_POS, 
               HIGH_INF_POS = CSQ_info$HIGH_INF_POS,
               MOTIF_SCORE_CHANGE = CSQ_info$MOTIF_SCORE_CHANGE, 
               IMPACT = CSQ_info$IMPACT, PICK = CSQ_info$PICK, 
               VARIANT_CLASS = CSQ_info$VARIANT_CLASS, TSL = CSQ_info$TSL, 
               HGVS_OFFSET = CSQ_info$HGVS_OFFSET, PHENO = CSQ_info$PHENO, 
               MINIMISED = '.', gnomAD_AF = CSQ_info$gnomAD_AF,
               gnomAD_AFR_AF = CSQ_info$gnomAD_AFR_AF, 
               gnomAD_AMR_AF = CSQ_info$gnomAD_AMR_AF,
               gnomAD_EAS_AF = CSQ_info$gnomAD_EAS_AF, 
               gnomAD_FIN_AF = CSQ_info$gnomAD_FIN_AF, 
               gnomAD_NFE_AF = CSQ_info$gnomAD_NFE_AF, 
               gnomAD_OTH_AF = CSQ_info$gnomAD_OTH_AF, 
               gnomAD_SAS_AF = CSQ_info$gnomAD_SAS_AF, 
               GENE_PHENO = CSQ_info$GENE_PHENO, FILTER = vcf_main$FILTER, 
               VAF = maf[, 47])
  maf <- cbind(maf, vcf_additional)
  
  ## add column isIndelAround
  maf <- cbind(maf, isIndelAround = 0)
  maf$Start_Position <- as.numeric(maf$Start_Position)
  maf$End_Position <- as.numeric(maf$End_Position)
  
  maf_arr <- arrange(maf, Chromosome, Variant_Type, Start_Position)
  chroms <- unique(maf_arr$Chromosome)
  for (c in 1:length(chroms)) {
    mafdat <- maf_arr[which(maf_arr$Chromosome == chroms[c]), ]
    if(any(mafdat$Variant_Type %in% c('INS', 'DEL'))) {
      mafIndel <- mafdat[which(mafdat$Variant_Type %in% c('INS', 'DEL')), ]
      mafSNP <- mafdat[which(!(mafdat$Variant_Type %in% c('INS', 'DEL'))), ]
      if(nrow(mafSNP) != 0) {
        for (i in 1:nrow(mafSNP)){
          if(any(abs(mafSNP$Start_Position[i] - c(mafIndel$End_Position, 
                                                  mafIndel$Start_Position)) <= 5) |
             any(abs(mafSNP$End_Position[i] - c(mafIndel$End_Position, 
                                                mafIndel$Start_Position)) <= 5)) {
            maf_arr[which(rownames(maf_arr) == rownames(mafSNP[i, ])), 
                    'isIndelAround'] <- 1
          }
        }
      }
      
    }
  }  
  
  
  if (writeFile) {
    message('The generated MAF file has been saved.')
    write.table(maf_arr, paste0(MAFdir, MAFfile), sep = "\t", 
                quote = FALSE, row.names = FALSE)
    
  } else{
    return(maf_arr)
  }
  
}





