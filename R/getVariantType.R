## get variant type and inframe info
getVariantType <- function(vcf_pos, ref, alt) {
  
  vcf_pos <- as.numeric(vcf_pos)
  # result that will contain: endPos, variantType
  if (nchar(ref) == nchar(alt) & (ref != '-') & (alt != '-')) {
    if (nchar(ref) == 1) {
      return(c(vcf_pos, 'SNP', TRUE))
    } else if (nchar(ref) == 2) {
      return(c(vcf_pos + 1, 'DNP', TRUE))
    } else if (nchar(ref) == 3) {
      return(c(vcf_pos + 2, 'TNP', TRUE))
    } else {
      return(c(vcf_pos + nchar(ref) - 1, 'ONP', TRUE))
    }
  } else if (nchar(ref) < nchar(alt) & (ref == '-')) {
    return(c(vcf_pos + nchar(alt) - 1, 'INS', (nchar(alt) %% 3 == 0)))
  } else if (nchar(ref) < nchar(alt) & (ref != '-')){
    return(c(vcf_pos + nchar(alt) - 1, 'INS', 
             ((nchar(alt) - nchar(ref)) %% 3 == 0)))
  } else if (nchar(ref) > nchar(alt) & (ref == '-')) {
    return(c(vcf_pos + nchar(ref) - 1, 'DEL', (nchar(ref) %% 3 == 0)))
  } else if (nchar(ref) > nchar(alt) & (ref != '-')){
    return(c(vcf_pos + nchar(ref) - 1, 'DEL', 
             ((nchar(ref) - nchar(alt)) %% 3 == 0)))
  }
}
