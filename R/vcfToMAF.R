vcfToMAF <- function(vcfFile, mAFfile = 'MAF.maf',MAFdir = './', 
                     tumorID = 'Extracted', normalID = 'Extracted', 
                     ncbiBuild = 'GRCh37', MAFcenter = '.', 
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

  
  
}





