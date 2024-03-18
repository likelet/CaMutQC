## select the proper transcript based on biotype, consequence and 
## length of transcript
selectMut <- function(charMatrix) { 
    if (nrow(charMatrix) == 1 ) { return(1)
    } else {
        numMatrix <- as.data.frame(matrix(ncol=1, nrow=nrow(charMatrix)))
        numMatrix[, 1] <- vapply(X=charMatrix[, 1], FUN=helperBiotype, FUN.VALUE=0)
        Biotypefreq <- data.frame(table(numMatrix[, 1]))
        Biotypefreq <- Biotypefreq[order(Biotypefreq$Var1, decreasing=FALSE), ]
        if (Biotypefreq[1, 2] == 1){
          # return the one with the highest biofunction priority
          return(as.numeric(rownames(charMatrix[numMatrix$V1 == Biotypefreq[1, 1], ])))
        } else {
          # keep rows with the highest biofunction priority, sort by consequence
          remMatrix <- charMatrix[numMatrix[, 1] == Biotypefreq[1, 1], ]
          remNumMatrix <- as.data.frame(matrix(ncol=1, nrow=nrow(remMatrix)))
          rownames(remNumMatrix) <- rownames(remMatrix)
          remNumMatrix[, 1] <- vapply(X=remMatrix[, 2], FUN=helperConsequence, FUN.VALUE=0)
          conseqFreq <- data.frame(table(remNumMatrix[, 1]))
          conseqFreq <- conseqFreq[order(conseqFreq$Var1, decreasing=FALSE), ]
        if (conseqFreq[1, 2] == 1){
            # return the one with the highest consequence priority
            return(as.numeric(rownames(remMatrix[remNumMatrix[, 1] == conseqFreq[1, 1], ])))
        } else {
            # work on length
            finalMatrix <- remMatrix[remNumMatrix[, 1] == conseqFreq[1, 1], ]
            ## choose the first transcript if the cDNA info of all variants are the same
            if (length(unique(finalMatrix$cDNA_position)) == 1) {
                return(as.numeric(rownames(finalMatrix)[1]))
            ## choose the only one transcript with cDNA position info
            } else if (sum(finalMatrix$cDNA_position != 'Missing') == 1) {
                return(as.numeric(rownames(finalMatrix[finalMatrix[, 3] != 'Missing', ])))
            } else {
            ## select the longest transcript
            finalMatrix <- finalMatrix[finalMatrix$cDNA_position != 'Missing', ]
            finalNumMatrix <- as.data.frame(matrix(ncol=1, nrow=nrow(finalMatrix)))
            rownames(finalNumMatrix) <- rownames(finalMatrix)
            # if "/" or "-" exists in the cDNA column, process this column
            # check "/" first
            if (any(str_detect(finalMatrix$cDNA_position, pattern="/"))){
                finalNumMatrix[, 1] <- as.numeric(str_split_fixed(finalMatrix[, 3], "/", 2)[,2])
                selectedTranscript <- as.numeric(rownames(
                finalMatrix[which.max(finalNumMatrix$V1), ]))
            # look for "-" if there is no "/"
            }else if (any(str_detect(finalMatrix$cDNA_position, pattern="-"))) {
                finalNumMatrix[, 1] <- as.numeric(str_split_fixed(finalMatrix[, 3], "-", 2)[,2])
                selectedTranscript <- as.numeric(rownames(
                finalMatrix[which.max(finalNumMatrix$V1), ]))
            }else{
                selectedTranscript <- as.numeric(rownames(
                finalMatrix[which.max(finalMatrix$cDNA_position), ]))
            }
            return(selectedTranscript)
        }
      }
    }
  }
}

## order biotype
## reference: vcf2maf, https://www.gencodegenes.org/pages/biotypes.html
getBiotypePriority <- function(biotype) {
    switch(biotype,
           'protein_coding' = 1, 'LRG_gene' = 2,
           'IG_C_gene' = 2, 'IG_D_gene' = 2, 'IG_J_gene' = 2, 'IG_LV_gene' = 2, 
           'IG_V_gene' = 2, 'TR_C_gene' = 2, 'TR_D_gene' = 2, 'TR_J_gene' = 2, 
           'TR_V_gene' = 2, 'miRNA' = 3, 'snRNA' = 3, 'snoRNA' = 3, 
           'ribozyme' = 3, 'tRNA' = 3, 'sRNA' = 3, 'scaRNA' = 3, 'scRNA' = 3, 
           'rRNA' = 3, 'lincRNA' = 3, 'lncRNA' = 3,
           'bidirectional_promoter_lncrna' = 3, 
           'bidirectional_promoter_lncRNA' = 3, 'known_ncrna' = 4,
           'vaultRNA' = 4,  'macro_lncRNA' = 4, 'Mt_tRNA' = 4, 'Mt_rRNA' = 4, 
           'antisense' = 5, 'antisense_RNA' = 5,
           'sense_intronic' = 5, 'sense_overlapping' = 5, 
           '3prime_overlapping_ncrna' = 5,  '3prime_overlapping_ncRNA' = 5, 
           'misc_RNA' = 5, 'vault_RNA' = 5, 'non_coding' = 5, 
           'regulatory_region' = 6, 'processed_transcript' = 6, 
           'protein_coding_CDS_not_defined' = 6, 'TEC' = 6, 
           'TF_binding_site' = 7, 'CTCF_binding_site' = 7, 
           'promoter_flanking_region' = 7,'enhancer' = 7, 'promoter' = 7, 
           'open_chromatin_region' = 7, 'retained_intron' = 7, 
           'nonsense_mediated_decay' = 7, 'non_stop_decay' = 7, 
           'ambiguous_orf' = 7, 'pseudogene' = 8, 'processed_pseudogene' = 8, 
           'polymorphic_pseudogene' = 8, 'protein_coding_LoF' = 8, 
           'retrotransposed' = 8, 'translated_processed_pseudogene' = 8, 
           'translated_unprocessed_pseudogene' = 8, 
           'transcribed_processed_pseudogene' = 8, 
           'transcribed_unprocessed_pseudogene' = 8, 
           'transcribed_unitary_pseudogene' = 8, 
           'transcribed_pseudogene' = 8, 
           'unitary_pseudogene' = 8, 'unprocessed_pseudogene' = 8, 
           'Mt_tRNA_pseudogene' = 8, 'tRNA_pseudogene' = 8, 
           'snoRNA_pseudogene' = 8, 'snRNA_pseudogene' = 8, 
           'scRNA_pseudogene' = 8,  'rRNA_pseudogene' = 8, 
           'misc_RNA_pseudogene' = 8, 'miRNA_pseudogene' = 8, 
           'IG_C_pseudogene' = 8, 'IG_D_pseudogene' = 8, 
           'IG_J_pseudogene' = 8, 'IG_V_pseudogene' = 8, 
           'TR_J_pseudogene' = 8, 'TR_V_pseudogene' = 8, 
           'artifact' = 9, 'Missing' = 10, 
    ) 
}

## order consequence
## reference: vcf2maf
## https://ensembl.org/info/genome/variation/prediction/predicted_data.html
getConsequencePriority <- function(consequence) {
    switch(consequence,
         'transcript_ablation' = 1, 'exon_loss_variant' = 1,
         'splice_donor_variant' = 2, 'splice_acceptor_variant' = 2,
         'stop_gained' = 3, 'frameshift_variant' = 3, 'stop_lost' = 3,
         'start_lost' = 4, 'initiator_codon_variant' = 4,
         'disruptive_inframe_insertion' = 5, 'disruptive_inframe_deletion' = 5,
         'inframe_insertion' = 5, 'inframe_deletion' = 5,
         'protein_altering_variant' = 5, 'missense_variant' = 6,
         'conservative_missense_variant' = 6, 'rare_amino_acid_variant' = 6,
         'transcript_amplification' = 7, 'splice_region_variant' = 8,
         'splice_donor_5th_base_variant' = 8, 
         'splice_donor_region_variant' = 8, 
         'splice_polypyrimidine_tract_variant' = 8, 
         'stop_retained_variant' = 9, 'start_retained_variant' = 9, 
         'synonymous_variant' = 9, 'incomplete_terminal_codon_variant' = 10,
         'coding_sequence_variant' = 11, 'mature_miRNA_variant' = 11,
         'exon_variant' = 11, '5_prime_UTR_variant' = 12,
         '5_prime_UTR_premature_start_codon_gain_variant' = 12,
         '3_prime_UTR_variant' = 12, 'non_coding_exon_variant' = 13,
         'non_coding_transcript_exon_variant' = 13,
         'non_coding_transcript_variant' = 14,
         'nc_transcript_variant' = 14, 'intron_variant' = 14,
         'intragenic_variant' = 14, 'INTRAGENIC' = 14, 
         'NMD_transcript_variant' = 15, 'upstream_gene_variant' = 16, 
         'downstream_gene_variant' = 16, 'TFBS_ablation' = 17,
         'TFBS_amplification' = 17,
         'TF_binding_site_variant' = 17,
         'regulatory_region_ablation' = 17,
         'regulatory_region_amplification' = 17,
         'regulatory_region_variant' = 17,
         'regulatory_region' = 17,
         'feature_elongation' = 18,
         'feature_truncation' = 18,
         'intergenic_variant' = 19,
         'intergenic_region' = 19,
         'Missing' = 20
    )
}

# helper function for getting the order of biotypes
helperBiotype <- function(biotype) {
    # handle situation when biotype is null
    if (biotype == '') {
        return(10)
    }else{
        # handle biotypes that are not included
        if (is.null(getBiotypePriority(biotype))) {
            return(10)
            warning("At least one biotype can not be recognized!")
        }else{
            return(getBiotypePriority(biotype))
        }
    }
}

# helper function for getting the order of consequences
helperConsequence <- function(consequence) {
    consequences <- strsplit(consequence, split="&")[[1]]
    conseqPriority <- 30
    # get the highest rank of consequences
    return(min(conseqPriority, vapply(consequences, FUN=getConsequencePriority, FUN.VALUE=1)))
}
