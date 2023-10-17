## select the proper transcript
selectMut <- function(charMatrix) {

  if (nrow(charMatrix) == 1 ) {
    return(1)
  } else {
    numMatrix <- as.data.frame(matrix(ncol = 2, nrow = nrow(charMatrix)))
    for (i in 1:nrow(numMatrix)) {
      numMatrix[i, 1] <- GetBiotypePriority(charMatrix[i, 1])
    }
    Biotypefreq <- data.frame(table(numMatrix[, 1]))
    Biotypefreq <- Biotypefreq[order(Biotypefreq$Var1, decreasing = FALSE), ]

    if (Biotypefreq[1, 2] == 1){

      # return the one with the highest biofunction priority
      return(as.numeric(rownames(numMatrix[which(numMatrix$V1
                                                 == Biotypefreq[1, 1]), ])))
    } else {

      # keep the rows with the highest biofunction priority, sort by consequence
      remMatrix <- charMatrix[which(numMatrix[, 1] == Biotypefreq[1, 1]), ]
      remNumMatrix <- as.data.frame(matrix(ncol = 2, nrow = nrow(remMatrix)))
      rownames(remNumMatrix) <- rownames(remMatrix)
      for (r in 1:nrow(remMatrix)) {
        Consequence <- strsplit(remMatrix[r, 2], split = "&")[[1]]
        conseqPriority <- 30
        for (c in Consequence) {
          conseqPriority <- min(conseqPriority, GetConsequencePriority(c))
        }
        remNumMatrix[r, 1] <- conseqPriority
      }

      Conseqfreq <- data.frame(table(remNumMatrix[, 1]))
      Conseqfreq <- Conseqfreq[order(Conseqfreq$Var1, decreasing = FALSE), ]

      if (Conseqfreq[1, 2] == 1){

        # return the one with the highest consequence priority
        return(as.numeric(rownames(remNumMatrix[which(remNumMatrix[, 1]
                                          == Conseqfreq[1, 1]), ])))
      } else {

        # work on length
        finalMatrix <- remMatrix[which(remNumMatrix[, 1]
                                       == Conseqfreq[1, 1]), ]
        #rownames(finalMatrix) <- rownames(remMatrix)

        ## choose the first transcript if all cDNA position info is missing
        if (all(finalMatrix$cDNA_position == 'Missing') |
            length(unique(finalMatrix$cDNA_position)) == 1) {
          return(as.numeric(rownames(finalMatrix)[1]))

          ## choose the only one transcript with cDNA position info
        } else if (sum(finalMatrix$cDNA_position != 'Missing') == 1) {
          return(as.numeric(rownames(finalMatrix[which(finalMatrix[, 3]
                                                       != 'Missing'), ])))
        } else {

          ## length of transcript
          finalMatrix <- finalMatrix[which(finalMatrix$cDNA_position
                                           != 'Missing'), ]

          if (length(grep('/', finalMatrix$cDNA_position))){
            Toplength <- as.numeric(strsplit(finalMatrix[1, 3], "/")[[1]][2])
            selectedTranscript <- as.numeric(rownames(finalMatrix[1, ]))
            #concurrent <- FALSE
            for (l in 2:nrow(finalMatrix)) {
              Tend <- as.numeric(strsplit(finalMatrix[l, 3], "/")[[1]][2])
              if (Tend > Toplength) {
                Toplength <- Tend
                selectedTranscript <- as.numeric(rownames(finalMatrix[l, ]))
                #concurrent <- FALSE
              }
            }
          }else if (length(grep('-', finalMatrix$cDNA_position))){
            Toplength <- as.numeric(strsplit(finalMatrix[1, 3], "-")[[1]][2])
            selectedTranscript <- as.numeric(rownames(finalMatrix[1, ]))
            #concurrent <- FALSE
            for (l in 2:nrow(finalMatrix)) {
              Tend <- as.numeric(strsplit(finalMatrix[l, 3], "-")[[1]][2])
              if (Tend > Toplength) {
                Toplength <- Tend
                selectedTranscript <- as.numeric(rownames(finalMatrix[l, ]))
                #concurrent <- FALSE
              }
            }
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
GetBiotypePriority <- function(biotype) {
  switch(biotype,
         'protein_coding' = 1,
         # Contains an open reading frame (ORF)
         'LRG_gene' = 2,
         # Gene in a "Locus Reference Genomic" region known to have
         # disease-related sequence variations
         'IG_C_gene' = 2,
         # Immunoglobulin (Ig) variable chain genes imported or
         # annotated according to the IMGT
         'IG_D_gene' = 2,
         # Immunoglobulin (Ig) variable chain genes imported or
         # annotated according to the IMGT
         'IG_J_gene' = 2,
         # Immunoglobulin (Ig) variable chain genes imported or
         # annotated according to the IMGT
         'IG_LV_gene' = 2, # Immunoglobulin (Ig) variable chain genes imported
         # or annotated according to the IMGT
         'IG_V_gene' = 2, # Immunoglobulin (Ig) variable chain genes imported
         # or annotated according to the IMGT
         'TR_C_gene' = 2, # T-cell receptor (TcR) genes imported or annotated
         # according to the IMGT
         'TR_D_gene' = 2, # T-cell receptor (TcR) genes imported or annotated
         # according to the IMGT
         'TR_J_gene' = 2, # T-cell receptor (TcR) genes imported or annotated
         # according to the IMGT
         'TR_V_gene' = 2, # T-cell receptor (TcR) genes imported or annotated
         # according to the IMGT
         'miRNA' = 3, # Non-coding RNA predicted using sequences
         # from RFAM and miRBase
         'snRNA' = 3, # Non-coding RNA predicted using sequences
         # from RFAM and miRBase
         'snoRNA' = 3, # Non-coding RNA predicted using sequences
         # from RFAM and miRBase
         'ribozyme' = 3, # Non-coding RNA predicted using sequences
         # from RFAM and miRBase
         'tRNA' = 3, #Added by Y. Boursin
         'sRNA' = 3, # Non-coding RNA predicted using sequences
         # from RFAM and miRBase
         'scaRNA' = 3, # Non-coding RNA predicted using sequences
         # from RFAM and miRBase
         'rRNA' = 3, # Non-coding RNA predicted using sequences
         # from RFAM and miRBase
         'lincRNA' = 3, # Long, intervening noncoding (linc) RNAs,
         # that can be found in evolutionarily conserved, intergenic regions
         'lncRNA' = 3,
         'bidirectional_promoter_lncrna' = 3, # A non-coding locus
         # that originates from within the promoter region of
         # a protein-coding gene, with transcription proceeding in the
         # opposite direction on the other strand
         'bidirectional_promoter_lncRNA' = 3, # A non-coding locus that
         # originates from within the promoter region of a protein-coding gene,
         # with transcription proceeding in the opposite direction
         # on the other strand
         'known_ncrna' = 4,
         'vaultRNA' = 4, # Short non coding RNA genes that form part
         # of the vault ribonucleoprotein complex
         'macro_lncRNA' = 4, # unspliced lncRNAs that are several kb in size
         'Mt_tRNA' = 4, # Non-coding RNA predicted using sequences
         # from RFAM and miRBase
         'Mt_rRNA' = 4, # Non-coding RNA predicted using sequences
         # from RFAM and miRBase
         'antisense' = 5, # Has transcripts that overlap the genomic span
         'lncRNA' = 5, # Generic long non-coding RNA biotype that replaced the following biotypes: 
         # 3prime_overlapping_ncRNA, antisense, bidirectional_promoter_lncRNA, 
         # lincRNA, macro_lncRNA, non_coding, processed_transcript, 
         # sense_intronic and sense_overlapping.
         # (exon or introns) of a protein-coding locus on the opposite strand
         'antisense_RNA' = 5, # Alias for antisense (Y. Boursin)
         'sense_intronic' = 5, # Long non-coding transcript in introns
         # of a coding gene that does not overlap any exons
         'sense_overlapping' = 5, # Long non-coding transcript that
         # contains a coding gene in its intron on the same strand
         '3prime_overlapping_ncrna' = 5, # Transcripts where ditag
         # and/or published experimental data strongly supports the existence
         # of short non-coding transcripts transcribed from the 3'UTR
         '3prime_overlapping_ncRNA' = 5, # Transcripts where ditag and/or published experimental data strongly supports the existence of short non-coding transcripts transcribed from the 3'UTR
         'misc_RNA' = 5, # Non-coding RNA predicted using sequences from RFAM and miRBase,
         'vault_RNA' = 5, # Short non coding RNA genes that form part of the vault ribonucleoprotein complex.
         'non_coding' = 5, # Transcript which is known from the literature to not be protein coding
         'regulatory_region' = 6, # A region of sequence that is involved in the control of a biological process
         'disrupted_domain' = 6, # Otherwise viable coding region omitted from this alternatively spliced transcript because the splice variation affects a region coding for a protein domain
         'processed_transcript' = 6, # Doesn't contain an ORF
         'protein_coding_CDS_not_defined' = 6, # Transcript that belongs to a protein_coding gene and doesn't contain an ORF. Replaces the processed_transcript transcript biotype in protein_coding genes.
         'TEC' = 6, # To be Experimentally Confirmed. This is used for non-spliced EST clusters that have polyA features. This category has been specifically created for the ENCODE project to highlight regions that could indicate the presence of protein coding genes that require experimental validation, either by 5' RACE or RT-PCR to extend the transcripts, or by confirming expression of the putatively-encoded peptide with specific antibodies
         'TF_binding_site' = 7, # A region of a nucleotide molecule that binds a Transcription Factor or Transcription Factor complex
         'CTCF_binding_site' = 7, # A transcription factor binding site with consensus sequence CCGCGNGGNGGCAG, bound by CCCTF-binding factor
         'promoter_flanking_region' = 7, # A region immediately adjacent to a promoter which may or may not contain transcription factor binding sites
         'enhancer' = 7, # A cis-acting sequence that increases the utilization of (some) eukaryotic promoters, and can function in either orientation and in any location (upstream or downstream) relative to the promoter
         'promoter' = 7, # A regulatory_region composed of the TSS(s) and binding sites for TF_complexes of the basal transcription machinery
         'open_chromatin_region' = 7, # A DNA sequence that in the normal state of the chromosome corresponds to an unfolded, un-complexed stretch of double-stranded DNA
         'retained_intron' = 7, # Alternatively spliced transcript believed to contain intronic sequence relative to other, coding, variants
         'nonsense_mediated_decay' = 7, # If the coding sequence (following the appropriate reference) of a transcript finishes >50bp from a downstream splice site then it is tagged as NMD. If the variant does not cover the full reference coding sequence then it is annotated as NMD if NMD is unavoidable i.e. no matter what the exon structure of the missing portion is the transcript will be subject to NMD
         'non_stop_decay' = 7, # Transcripts that have polyA features (including signal) without a prior stop codon in the CDS, i.e. a non-genomic polyA tail attached directly to the CDS without 3' UTR. These transcripts are subject to degradation
         'ambiguous_orf' = 7, # Transcript believed to be protein coding, but with more than one possible open reading frame
         'pseudogene' = 8, # Have homology to proteins but generally suffer from a disrupted coding sequence and an active homologous gene can be found at another locus. Sometimes these entries have an intact coding sequence or an open but truncated ORF, in which case there is other evidence used (for example genomic polyA stretches at the 3' end) to classify them as a pseudogene. Can be further classified as one of the following
         'processed_pseudogene' = 8, # Pseudogene that lack introns and is thought to arise from reverse transcription of mRNA followed by reinsertion of DNA into the genome
         'polymorphic_pseudogene' = 8, # Pseudogene owing to a SNP/DIP but in other individuals/haplotypes/strains the gene is translated
         'protein_coding_LoF' = 8, # Not translated in the reference genome owing to a SNP/DIP but in other individuals/haplotypes/strains the transcript is translated. Replaces the polymorphic_pseudogene transcript biotype.
         'retrotransposed' = 8, # Pseudogene owing to a reverse transcribed and re-inserted sequence
         'translated_processed_pseudogene' = 8, # Pseudogenes that have mass spec data suggesting that they are also translated
         'translated_unprocessed_pseudogene' = 8, # Pseudogenes that have mass spec data suggesting that they are also translated
         'transcribed_processed_pseudogene' = 8, # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
         'transcribed_unprocessed_pseudogene' = 8, # Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
         'transcribed_unitary_pseudogene' = 8, #Pseudogene where protein homology or genomic structure indicates a pseudogene, but the presence of locus-specific transcripts indicates expression
         'unitary_pseudogene' = 8, # A species specific unprocessed pseudogene without a parent gene, as it has an active orthologue in another species
         'unprocessed_pseudogene' = 8, # Pseudogene that can contain introns since produced by gene duplication
         'Mt_tRNA_pseudogene' = 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
         'tRNA_pseudogene' = 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
         'snoRNA_pseudogene' = 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
         'snRNA_pseudogene' = 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
         'scRNA_pseudogene' = 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
         'rRNA_pseudogene' = 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
         'misc_RNA_pseudogene' = 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
         'miRNA_pseudogene' = 8, # Non-coding RNAs predicted to be pseudogenes by the Ensembl pipeline
         'IG_C_pseudogene' = 8, # Inactivated immunoglobulin gene
         'IG_D_pseudogene' = 8, # Inactivated immunoglobulin gene
         'IG_J_pseudogene' = 8, # Inactivated immunoglobulin gene
         'IG_V_pseudogene' = 8, # Inactivated immunoglobulin gene
         'TR_J_pseudogene' = 8, # Inactivated immunoglobulin gene
         'TR_V_pseudogene' = 8, # Inactivated immunoglobulin gene
         'artifact' = 9, # Used to tag mistakes in the public databases (Ensembl/SwissProt/Trembl)
         'Missing' = 10) # Used to process missing value
}


## order consequence
GetConsequencePriority <- function(consequence) {
  switch(consequence,
         'transcript_ablation' = 1,
         # A feature ablation whereby the deleted region includes a transcript feature
         'exon_loss_variant' = 1,
         # A sequence variant whereby an exon is lost from the transcript
         'splice_donor_variant' = 2,
         # A splice variant that changes the 2 base region at the 5' end of an intron
         'splice_acceptor_variant' = 2,
         # A splice variant that changes the 2 base region at the 3' end of an intron
         'stop_gained' = 3,
         # A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript
         'frameshift_variant' = 3,
         # A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three
         'stop_lost' = 3,
         # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
         'start_lost' = 4,
         # A codon variant that changes at least one base of the canonical start codon
         'initiator_codon_variant' = 4,
         # A codon variant that changes at least one base of the first codon of a transcript
         'disruptive_inframe_insertion' = 5,
         # An inframe increase in cds length that inserts one or more codons into the coding sequence within an existing codon
         'disruptive_inframe_deletion' = 5,
         # An inframe decrease in cds length that deletes bases from the coding sequence starting within an existing codon
         'inframe_insertion' = 5,
         # An inframe non synonymous variant that inserts bases into the coding sequence
         'inframe_deletion' = 5,
         # An inframe non synonymous variant that deletes bases from the coding sequence
         'protein_altering_variant' = 5,
         # A sequence variant which is predicted to change the protein encoded in the coding sequence
         'missense_variant' = 6,
         # A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved
         'conservative_missense_variant' = 6,
         # A sequence variant whereby at least one base of a codon is changed resulting in a codon that encodes for a different but similar amino acid. These variants may or may not be deleterious
         'rare_amino_acid_variant' = 6,
         # A sequence variant whereby at least one base of a codon encoding a rare amino acid is changed, resulting in a different encoded amino acid
         'transcript_amplification' = 7,
         # A feature amplification of a region containing a transcript
         'splice_region_variant' = 8,
         # A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron
         'stop_retained_variant' = 9,
         # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
         'synonymous_variant' = 9,
         # A sequence variant where there is no resulting change to the encoded amino acid
         'incomplete_terminal_codon_variant' = 10,
         # A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed
         'coding_sequence_variant' = 11,
         # A sequence variant that changes the coding sequence
         'mature_miRNA_variant' = 11,
         # A transcript variant located with the sequence of the mature miRNA
         'exon_variant' = 11,
         # A sequence variant that changes exon sequence
         '5_prime_UTR_variant' = 12, # A UTR variant of the 5' UTR
         '5_prime_UTR_premature_start_codon_gain_variant' = 12,
         # snpEff-specific effect, creating a start codon in 5' UTR
         '3_prime_UTR_variant' = 12, # A UTR variant of the 3' UTR
         'non_coding_exon_variant' = 13,
         # A sequence variant that changes non-coding exon sequence
         'non_coding_transcript_exon_variant' = 13,
         # snpEff-specific synonym for non_coding_exon_variant
         'non_coding_transcript_variant' = 14,
         # A transcript variant of a non coding RNA gene
         'nc_transcript_variant' = 14,
         # A transcript variant of a non coding RNA gene (older alias for non_coding_transcript_variant)
         'intron_variant' = 14,
         # A transcript variant occurring within an intron
         'intragenic_variant' = 14,
         # A variant that occurs within a gene but falls outside of all transcript features. This occurs when alternate transcripts of a gene do not share overlapping sequence
         'INTRAGENIC' = 14, # snpEff-specific synonym of intragenic_variant
         'NMD_transcript_variant' = 15,
         # A variant in a transcript that is the target of NMD
         'upstream_gene_variant' = 16, # A sequence variant located 5' of a gene
         'downstream_gene_variant' = 16,
         # A sequence variant located 3' of a gene
         'TFBS_ablation' = 17,
         # A feature ablation whereby the deleted region includes a transcription factor binding site
         'TFBS_amplification' = 17,
         # A feature amplification of a region containing a transcription factor binding site
         'TF_binding_site_variant' = 17,
         # A sequence variant located within a transcription factor binding site
         'regulatory_region_ablation' = 17,
         # A feature ablation whereby the deleted region includes a regulatory region
         'regulatory_region_amplification' = 17,
         # A feature amplification of a region containing a regulatory region
         'regulatory_region_variant' = 17,
         # A sequence variant located within a regulatory region
         'regulatory_region' = 17,
         # snpEff-specific effect that should really be regulatory_region_variant
         'feature_elongation' = 18,
         # A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence
         'feature_truncation' = 18,
         # A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence
         'intergenic_variant' = 19,
         # A sequence variant located in the intergenic region, between genes
         'intergenic_region' = 19,
         # snpEff-specific effect that should really be intergenic_variant
         'Missing' = 20
         # Used to process missing value
  )
}



