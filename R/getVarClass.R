getVarClass <- function(vcf_consequence, variantType, inframe) {
  if (vcf_consequence %in% c('splice_acceptor_variant', 'splice_donor_variant', 
                             'transcript_ablation', 'exon_loss_variant')) {
    return('Splice_Site')
  }else if (vcf_consequence == 'stop_gained') {return('Nonsense_Mutation')
  } else if (vcf_consequence == 'stop_lost') {return('Nonstop_Mutation')
  }else if(vcf_consequence == 'splice_region_variant') {return('Splice_Region')
  }else if(vcf_consequence == '3_prime_UTR_variant') {return('3\'UTR')
  }else if(vcf_consequence == 'upstream_gene_variant') {return('5\'Flank')
  }else if(vcf_consequence == 'downstream_gene_variant') {return('3\'Flank')
  }else if(vcf_consequence %in% c('transcript_amplification', 'INTRAGENIC',
                                  'intron_variant', 'intragenic_variant')) {
    return('Intron')
  }else if (vcf_consequence %in% c('incomplete_terminal_codon_variant', 
                                   'synonymous_variant','stop_retained_variant', 
                                   'NMD_transcript_variant')) {return('Silent')
  }else if(vcf_consequence %in% 
           c('mature_miRNA_variant','nc_transcript_variant', 'exon_variant', 
              'non_coding_exon_variant', 'non_coding_transcript_exon_variant', 
              'non_coding_transcript_variant')) {return('RNA')
  }else if(vcf_consequence %in% c('5_prime_UTR_variant', 
        '5_prime_UTR_premature_start_codon_gain_variant')) {return('5\'UTR')
  }else if(vcf_consequence %in% 
           c('TF_binding_site_variant', 'regulatory_region_variant', 
             'regulatory_region', 'intergenic_variant', 'intergenic_region')) {
    return('IGR')
  }else if(vcf_consequence %in% 
           c('missense_variant', 'coding_sequence_variant', 
             'conservative_missense_variant', 'rare_amino_acid_variant')) {
    return('Missense_Mutation')
  } else if (variantType == 'DEL' & (vcf_consequence == 'frameshift_variant' | 
            (vcf_consequence == 'protein_altering_variant' & (!inframe)))){
    return('Frame_Shift_Del')
  } else if (variantType == 'INS' & (vcf_consequence == 'frameshift_variant' | 
           (vcf_consequence == 'protein_altering_variant' & (!inframe)))) {
    return('Frame_Shift_Ins')
  } else if (vcf_consequence %in% c('initiator_codon_variant', 'start_lost')) {
    return('Translation_Start_Site')
  } else if ((vcf_consequence %in% c('inframe_insertion', 
                                     'disruptive_inframe_insertion')) | 
          (vcf_consequence == 'protein_altering_variant' & inframe 
          & variantType == 'INS')) {return('In_Frame_Ins')
  } else if ((vcf_consequence %in% c('inframe_deletion', 
                                     'disruptive_inframe_deletion')) | 
             (vcf_consequence == 'protein_altering_variant' & inframe 
              & variantType == 'DEL')) {return('In_Frame_Del')
  }else{return('Targeted_Region')}
} 