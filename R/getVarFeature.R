## get correct variant Position, Variant_Type, Ref allele and Alt allele
getVarFeature <- function(vcf_pos, ref, alt, csqalt) {
    vcf_pos <- as.numeric(vcf_pos)
    # first decide which alt is correct
    if (alt != csqalt){
    # if > one alt in exists, select one based on selected transcript (CSQ_info)
    # varscan: the alt field can have A/G, 
    # mutect and muse: the alt field may have A,G
        if (str_detect(alt, ",") | str_detect(alt, "/")){
          alts <- str_split(alt, "[,/]")[[1]]
          # handle special case: REF: C. Alt: T,CCAT. csqalt: CCAT
          if (csqalt %in% alts) {
            if (nchar(csqalt) == 1 & nchar(ref) == 1) { alt <- csqalt
            }else{
              ref <- remove1stString(ref)
              alt <- remove1stString(csqalt)}
          }else{
            # for a SNP, DNP, TNP or ONP, just select the first alt as final alt
            if (nchar(ref) == nchar(alts[1])) {
              if (csqalt %in% alts) {alt <- csqalt}else{alt <- alts[1]}
            }else{
              alts <- unlist(lapply(alts, remove1stString))
              ref <- remove1stString(ref)
              # if the alt in transcript can be found in alts, return it
              # otherwise, return the first alt in alts
              if (csqalt %in% alts) { alt <- csqalt } else {alt <- alts[1]} }
          }
        }else if (nchar(ref) != nchar(alt)){
          # remove the same prefix of ref and alt
          ref <- remove1stString(ref)
          alt <- remove1stString(alt) }
    }
    # SNV. If ref and alt have the same length and there is no -
    if (nchar(ref) == nchar(alt) & !("-" %in% c(ref, alt))){
        if (nchar(ref) == 1) {
          return(list(vcf_pos, vcf_pos, 'SNP', ref, alt, TRUE))
        } else if (nchar(ref) == 2) {
          return(list(vcf_pos, vcf_pos + 1, 'DNP', ref, alt, TRUE))
        } else if (nchar(ref) == 3) {
          return(list(vcf_pos, vcf_pos + 2, 'TNP', ref, alt, TRUE))
        } else {
          return(list(vcf_pos,vcf_pos + nchar(ref) - 1, 'ONP', ref, alt, TRUE))}
    } else if (ref == "-" & alt != "-"){
        # INS
        start_pos <- vcf_pos + 1
        return(list(start_pos, start_pos + 1, "INS", ref, alt, 
                    (nchar(alt)) %% 3 == 0))
    }else if (alt == "-" & ref != "-"){
        # DEL
        start_pos <- vcf_pos + 1
        return(list(start_pos, start_pos + nchar(ref) - 1, "DEL", ref, alt, 
                    (nchar(ref)) %% 3 == 0))
        # if the first bp is still the same, remove the first character
    }else if (substring(ref, 1, 1) == substring(alt, 1, 1)){
        ref <- remove1stString(ref)
        alt <- remove1stString(alt)
        return(list(vcf_pos + 2, vcf_pos + 2, "INS", ref, alt, 
                    (nchar(alt)) %% 3 == 0))
        # INS, handle cases: Ref: C. Alt: AC. csqalt: CAC
    }else if (nchar(alt) > nchar(ref)){
        start_pos <- vcf_pos + 1
        ref <- "-"
        alt <- substring(alt,1,1)
        return(list(start_pos, start_pos + 1, "INS", ref, alt, 
                    (nchar(alt)) %% 3 == 0))
    }else{
        mes <- paste0("Ref: ", ref, ". Alt: ", alt, ". csqalt: ", csqalt)
        stop(mes)
    }
}

# remove the first string 
remove1stString <- function(chars) {
  # move the 1st character
  remain_char <- substring(chars, 2, nchar(chars))
  # return - if no characters remain
  if (remain_char == "") {
    return("-")
  }else{
    return(remain_char)
  }
}

