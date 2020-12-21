library(vcfR)
library(org.Hs.eg.db)
library(clusterProfiler)

vcfToMAF <- function(vcfFile, mAFfile = 'MAF.maf',MAFdir = './', 
                     tumorID = 'Extracted', normalID = 'Extracted', 
                     ncbiBuild = 'GRCh37', MAFcenter = '.', MAFstrand = '+', 
                     species = 'homo_sapiens'){
  
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
                     'Chromosome', 'Start_Position', 'End_Position', 
                     'Strand', 'Variant_Classification', 'Variant_Type',
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
  CSQ_headers <- strsplit(INFOFrame[25, 4], 
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
  if (tumorID == 'Extracted'){
    tumorSampleLine <- vcf_header$`Anno_Vcf@meta`[grep('tumor_sample=', 
                                                       vcf_header$`Anno_Vcf@meta`)]
    tumorID <- strsplit(tumorSampleLine, split = '=')[[1]][2]
    
  }
  
  if (normalID == 'Extracted'){
    normalSampleLine <- vcf_header$`Anno_Vcf@meta`[grep('normal_sample=', 
                                                       vcf_header$`Anno_Vcf@meta`)]
    normalID <- strsplit(normalSampleLine, split = '=')[[1]][2]
  }
  
  ### check the consistence of tumorID/normalID with colnames of FORMAT column
  if (!(tumorID %in% colnames(vcf_additional))){
    stop('Tumor ID is inconsistent. 
     Set to \'Extracted\' if you want tumor ID to be extracted from VCF file')
  }else if(!(normalID %in% colnames(vcf_additional))){
    stop('Normal ID is inconsistent.
     Set to \'Extracted\' if you want normal ID to be extracted from VCF file')
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
  
  for (i in 1:nrow(maf)) {
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
    AD <- strsplit(strsplit(vcf_additional[i, tumorID], ":")[[1]][AD_loc], 
                   ",")[[1]][2]
    DP <- strsplit(vcf_additional[i, tumorID], ":")[[1]][DP_loc]
    maf[i, 'VAF'] <- as.numeric(AD)/as.numeric(DP)
  }
  
}





