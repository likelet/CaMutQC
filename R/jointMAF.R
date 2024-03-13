jointMAF <- function(mafDat, csqInfo, vcfMain){
    mafFinal <- cbind(mafDat, Allele = csqInfo$Allele, Gene = csqInfo$Gene, 
        Feature = csqInfo$Feature, Feature_type = csqInfo$Feature_type,
        One_Consequence = csqInfo$Consequence, Consequence = csqInfo$Consequence,
        cDNA_position = getValue(csqInfo$cDNA_position),
        CDS_position = getValue(csqInfo$CDS_position),
        Protein_position = getValue(csqInfo$Protein_position),
        Amino_acids = getValue(csqInfo$Amino_acids), 
        Codons = getValue(csqInfo$Codons),
        Existing_variation = getValue(csqInfo$Existing_variation),
        ALLELE_NUM = getValue(csqInfo$ALLELE_NUM),
        DISTANCE = getValue(csqInfo$DISTANCE), 
        TRANSCRIPT_STRAND = getValue(csqInfo$STRAND), 
        SYMBOL = getValue(csqInfo$SYMBOL), 
        SYMBOL_SOURCE = getValue(csqInfo$SYMBOL_SOURCE), 
        HGNC_ID = getValue(csqInfo$HGNC_ID), 
        BIOTYPE = getValue(csqInfo$BIOTYPE), 
        CANONICAL = getValue(csqInfo$CANONICAL),
        CCDS = getValue(csqInfo$CCDS), ENSP = getValue(csqInfo$ENSP), 
        SWISSPROT = getValue(csqInfo$SWISSPROT), TREMBL = getValue(csqInfo$TREMBL), 
        UNIPARC = getValue(csqInfo$UNIPARC), RefSeq = getValue(csqInfo$RefSeq), 
        SIFT = getValue(csqInfo$SIFT), PolyPhen = getValue(csqInfo$PolyPhen), 
        Exon = getValue(csqInfo$EXON), INTRON = getValue(csqInfo$INTRON), 
        DOMAINS = getValue(csqInfo$DOMAINS), GMAF = getValue(csqInfo$AF), 
        AFR_MAF = getValue(csqInfo$AFR_AF), AMR_MAF = getValue(csqInfo$AMR_AF), 
        ASN_MAF = getValue(csqInfo$ASN_MAF), EAS_MAF = getValue(csqInfo$EAS_AF),
        EUR_MAF = getValue(csqInfo$EUR_AF), SAS_MAF = getValue(csqInfo$SAS_AF), 
        AA_MAF = getValue(csqInfo$AA_AF), EA_MAF = getValue(csqInfo$EA_AF), 
        CLIN_SIG = getValue(csqInfo$CLIN_SIG), SOMATIC = getValue(csqInfo$SOMATIC), 
        PUBMED = getValue(csqInfo$PUBMED),
        MOTIF_NAME = getValue(csqInfo$MOTIF_NAME),
        MOTIF_POS = getValue(csqInfo$MOTIF_POS),
        HIGH_INF_POS = getValue(csqInfo$HIGH_INF_POS),
        MOTIF_SCORE_CHANGE = getValue(csqInfo$MOTIF_SCORE_CHANGE),
        IMPACT = getValue(csqInfo$IMPACT), PICK = getValue(csqInfo$PICK),
        VARIANT_CLAS = getValue(csqInfo$VARIANT_CLASS),
        TSL = getValue(csqInfo$TSL), HGVS_OFFSET = getValue(csqInfo$HGVS_OFFSET),
        PHENO = getValue(csqInfo$PHENO), MINIMISED = '', 
        gnomAD_AF = getGnominfo(csqInfo, fieldName="AF"),
        gnomAD_AFR_AF = getGnominfo(csqInfo, fieldName="AFR_AF"),
        gnomAD_AMR_AF = getGnominfo(csqInfo, fieldName="AMR_AF"),
        gnomAD_ASJ_AF = getGnominfo(csqInfo, fieldName="ASJ_AF"),
        gnomAD_EAS_AF = getGnominfo(csqInfo, fieldName="EAS_AF"),
        gnomAD_FIN_AF = getGnominfo(csqInfo, fieldName="FIN_AF"),
        gnomAD_NFE_AF = getGnominfo(csqInfo, fieldName="NFE_AF"),
        gnomAD_OTH_AF = getGnominfo(csqInfo, fieldName="OTH_AF"),
        gnomAD_SAS_AF = getGnominfo(csqInfo, fieldName="SAS_AF"),
        ExAC_AF = getValue(csqInfo$ExAC_AF),
        ExAC_AF_adj = getValue(csqInfo$ExAC_AF_adj),
        ExAC_AF_AFR = getValue(csqInfo$ExAC_AF_AFR),
        ExAC_AF_AMR = getValue(csqInfo$ExAC_AF_AMR),
        ExAC_AF_EAS = getValue(csqInfo$ExAC_AF_EAS),
        ExAC_AF_FIN = getValue(csqInfo$ExAC_AF_FIN),
        ExAC_AF_NFE = getValue(csqInfo$ExAC_AF_NFE),
        ExAC_AF_OTH = getValue(csqInfo$ExAC_AF_OTH),
        ExAC_AF_SAS = getValue(csqInfo$ExAC_AF_SAS),
        GENE_PHENO = getValue(csqInfo$GENE_PHENO), FILTER = vcfMain$FILTER)}

getValue <- function(chars){
    if (is.null(chars)){
        return("")
    }else{
        return(chars)
    }
}

# handle gnomADe and gnomAD cases
getGnominfo <- function(csqInfo, fieldName){
    # gnomAD changes to gnomADe in VEP v110 (at least)
    fieldName1 <- paste("gnomAD", fieldName, sep = "_")
    fieldName2 <- paste("gnomADe", fieldName, sep = "_")
    if (fieldName1 %in% colnames(csqInfo)){
        return(csqInfo[, fieldName1])
    }else if (fieldName2 %in% colnames(csqInfo)){
        return(csqInfo[, fieldName2])
    }else{
        return("")
    }
}

