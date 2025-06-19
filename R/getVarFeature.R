## get correct variant Position, Variant_Type, Ref allele and Alt allele
getVarFeature <- function(vcfPos, ref, alt, csqalt) {
    vcfPos <- as.numeric(vcfPos)
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
              ref <- remove1ststring(ref)
              alt <- remove1ststring(csqalt)}
          }else{
            # for a SNP, DNP, TNP or ONP, just select the first alt as final alt
            if (nchar(ref) == nchar(alts[1])) {
              if (csqalt %in% alts) {alt <- csqalt}else{alt <- alts[1]}
            }else{
              alts <- unlist(lapply(alts, remove1ststring))
              ref <- remove1ststring(ref)
              # if the alt in transcript can be found in alts, return it
              # otherwise, return the first alt in alts
              if (csqalt %in% alts) { alt <- csqalt } else {alt <- alts[1]} }
          }
        }else if (nchar(ref) != nchar(alt)){
          # remove the same prefix of ref and alt
          ref <- remove1ststring(ref)
          alt <- remove1ststring(alt) }
    }
    # SNV. If ref and alt have the same length and there is no -
    if (nchar(ref) == nchar(alt) & !("-" %in% c(ref, alt))){
        if (nchar(ref) == 1) {
          return(list(vcfPos, vcfPos, 'SNP', ref, alt, TRUE))
        } else if (nchar(ref) == 2) {
          return(list(vcfPos, vcfPos + 1, 'DNP', ref, alt, TRUE))
        } else if (nchar(ref) == 3) {
          return(list(vcfPos, vcfPos + 2, 'TNP', ref, alt, TRUE))
        } else {
          return(list(vcfPos, vcfPos + nchar(ref) - 1, 'ONP', ref, alt, TRUE))}
    } else if (ref == "-" & alt != "-"){
        # INS
        start_pos <- vcfPos + 1
        return(list(start_pos, start_pos + 1, "INS", ref, alt, 
                    (nchar(alt)) %% 3 == 0))
    }else if (alt == "-" & ref != "-"){
        # DEL
        start_pos <- vcfPos + 1
        return(list(start_pos, start_pos + nchar(ref) - 1, "DEL", ref, alt, 
                    (nchar(ref)) %% 3 == 0))
        # if the first bp is still the same, remove the first character
    }else if (substring(ref, 1, 1) == substring(alt, 1, 1)){
        ref <- remove1ststring(ref)
        alt <- remove1ststring(alt)
        return(list(vcfPos + 2, vcfPos + 2, "INS", ref, alt, 
                    (nchar(alt)) %% 3 == 0))
        # INS, handle cases: Ref: C. Alt: AC. csqalt: CAC
    }else if (nchar(alt) > nchar(ref)){
        start_pos <- vcfPos + 1
        ref <- "-"
        alt <- substring(alt, 1, 1)
        return(list(start_pos, start_pos + 1, "INS", ref, alt, 
                    (nchar(alt)) %% 3 == 0))
    # handle case: Ref: GTGG. Alt: TGG. csqalt: TGG
    }else if(nchar(alt) < nchar(ref) && any(grepl(alt, ref)) ) {
        ref_len <- nchar(ref)
        alt_len <- nchar(alt)
        # iterate through every possible case
        for (i in 1:(ref_len - alt_len + 1)) {
            if (substr(ref, i, i + alt_len - 1) == alt) {
                # calculate the deleted base
                del_start <- 1
                del_end <- i - 1
                deleted_bases <- substr(ref, del_start, del_end)
                length_changed <- del_end - del_start + 1
                # adjust ref and alt
                new_ref <- deleted_bases
                new_alt <- "-"
                # modify pos
                start_pos <- vcfPos
                end_pos <- vcfPos + length_changed - 1
                # return results
                return(list(start_pos, end_pos, 'DEL', new_ref, new_alt,
                            (length_changed %% 3 == 0)))
            }
        }
    }else{
        mes <- paste0("Ref: ", ref, ". Alt: ", alt, ". csqalt: ", csqalt)
        print("This is a case that CaMutQC cannot handle during MAF transformation.")
        print("Please contact the developer through email or post on github issue, thanks!")
        stop(mes)
    }
}

# remove the first string 
remove1ststring <- function(chars) {
    # move the 1st character
    remain_char <- substring(chars, 2, nchar(chars))
    # return - if no characters remain
    if (remain_char == "") {
        return("-")
    }else{
        return(remain_char)
    }
}

