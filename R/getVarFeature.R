## get correct variant Position, Variant_Type, Ref allele and Alt allele
getVarFeature <- function(vcf_pos, ref, alt) {
  vcf_pos <- as.numeric(vcf_pos)
  # SNV
  if (nchar(ref) == nchar(alt)){
    if (nchar(ref) == 1) {
      return(list(vcf_pos, vcf_pos, 'SNP', ref, alt, TRUE))
    } else if (nchar(ref) == 2) {
      return(list(vcf_pos, vcf_pos + 1, 'DNP', ref, alt, TRUE))
    } else if (nchar(ref) == 3) {
      return(list(vcf_pos, vcf_pos + 2, 'TNP', ref, alt, TRUE))
    } else {
      return(list(vcf_pos, vcf_pos + nchar(ref) - 1, 'ONP', ref, alt, TRUE))
    }
  } else if (nchar(ref) < nchar(alt)){
    # INS
    alt_allele <- substring(alt, nchar(ref)+1, nchar(alt))
    return(list(vcf_pos, vcf_pos + 1, "INS", "-",
                alt_allele, (nchar(alt) - nchar(ref)) %% 3 == 0))
  }else{
    # DEL
    start_pos <- vcf_pos + nchar(alt)
    ref_allele <- substring(ref, nchar(alt)+1, nchar(ref))
    return(list(start_pos,start_pos + nchar(ref) - nchar(alt) - 1, "DEL", 
                ref_allele, "-", (nchar(alt) - nchar(ref)) %% 3 == 0))
  }
}
