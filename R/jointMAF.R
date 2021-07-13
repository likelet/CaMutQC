jointMAF <- function(mafdf, CSQ_info, vcfMain){
  mafFinal <- cbind(mafdf, Allele = getValue(CSQ_info$Allele),
                    Gene = getValue(CSQ_info$Gene),
                    Feature = getValue(CSQ_info$Feature),
                    Feature_type = getValue(CSQ_info$Feature_type),
                    One_Consequence = getValue(CSQ_info$Consequence),
                    Consequence = getValue(CSQ_info$Consequence),
                    cDNA_position = getValue(CSQ_info$cDNA_position),
                    CDS_position = getValue(CSQ_info$CDS_position),
                    Protein_position = getValue(CSQ_info$Protein_position),
                    Amino_acids = getValue(CSQ_info$Amino_acids),
                    Codons = getValue(CSQ_info$Codons),
                    Existing_variation = getValue(CSQ_info$Existing_variation),
                    ALLELE_NUM = getValue(CSQ_info$ALLELE_NUM),
                    DISTANCE = getValue(CSQ_info$DISTANCE),
                    TRANSCRIPT_STRAND = getValue(CSQ_info$STRAND),
                    SYMBOL = getValue(CSQ_info$SYMBOL),
                    SYMBOL_SOURCE = getValue(CSQ_info$SYMBOL_SOURCE),
                    HGNC_ID = getValue(CSQ_info$HGNC_ID),
                    BIOTYPE = getValue(CSQ_info$BIOTYPE),
                    CANONICAL = getValue(CSQ_info$CANONICAL),
                    CCDS = getValue(CSQ_info$CCDS),
                    ENSP = getValue(CSQ_info$ENSP),
                    SWISSPROT = getValue(CSQ_info$SWISSPROT),
                    TREMBL = getValue(CSQ_info$TREMBL),
                    UNIPARC = getValue(CSQ_info$UNIPARC),
                    RefSeq = getValue(CSQ_info$RefSeq),
                    SIFT = getValue(CSQ_info$SIFT),
                    PolyPhen = getValue(CSQ_info$PolyPhen),
                    Exon = getValue(CSQ_info$EXON),
                    INTRON = getValue(CSQ_info$INTRON),
                    DOMAINS = getValue(CSQ_info$DOMAINS),
                    GMAF = getValue(CSQ_info$gnomAD_AF),
                    AFR_MAF = getValue(CSQ_info$AFR_AF),
                    AMR_MAF = getValue(CSQ_info$AMR_AF),
                    ASN_MAF = '.',
                    EAS_MAF = getValue(CSQ_info$EAS_AF),
                    EUR_MAF = getValue(CSQ_info$EUR_AF),
                    SAS_MAF = getValue(CSQ_info$SAS_AF),
                    AA_MAF = getValue(CSQ_info$AA_AF),
                    EA_MAF = getValue(CSQ_info$EA_AF),
                    CLIN_SIG = getValue(CSQ_info$CLIN_SIG),
                    SOMATIC = getValue(CSQ_info$SOMATIC),
                    PUBMED = getValue(CSQ_info$PUBMED),
                    MOTIF_NAME = getValue(CSQ_info$MOTIF_NAME),
                    MOTIF_POS = getValue(CSQ_info$MOTIF_POS),
                    HIGH_INF_POS = getValue(CSQ_info$HIGH_INF_POS),
                    MOTIF_SCORE_CHANGE = getValue(CSQ_info$MOTIF_SCORE_CHANGE),
                    IMPACT = getValue(CSQ_info$IMPACT),
                    PICK = getValue(CSQ_info$PICK),
                    VARIANT_CLASS = getValue(CSQ_info$VARIANT_CLASS),
                    TSL = getValue(CSQ_info$TSL),
                    HGVS_OFFSET = getValue(CSQ_info$HGVS_OFFSET),
                    PHENO = getValue(CSQ_info$PHENO),
                    MINIMISED = '.',
                    gnomAD_AF = getValue(CSQ_info$gnomAD_AF),
                    gnomAD_AFR_AF = getValue(CSQ_info$gnomAD_AFR_AF),
                    gnomAD_AMR_AF = getValue(CSQ_info$gnomAD_AMR_AF),
                    gnomAD_EAS_AF = getValue(CSQ_info$gnomAD_EAS_AF),
                    gnomAD_FIN_AF = getValue(CSQ_info$gnomAD_FIN_AF),
                    gnomAD_NFE_AF = getValue(CSQ_info$gnomAD_NFE_AF),
                    gnomAD_OTH_AF = getValue(CSQ_info$gnomAD_OTH_AF),
                    gnomAD_SAS_AF = getValue(CSQ_info$gnomAD_SAS_AF),
                    GENE_PHENO = getValue(CSQ_info$GENE_PHENO),
                    FILTER = vcfMain$FILTER)
  return(mafFinal)

}

getValue <- function(chars){
  if (is.null(chars)){
    return(NA)
  }else{
    return(chars)
  }
}