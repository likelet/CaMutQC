## PON filtering using external dataset and info flag

somatic38 <- read.vcfR('../PONdata/somatic-hg38_1000g_pon.hg38.vcf')
somatic37 <- read.vcfR('../PONdata/somatic-b37_Mutect2-exome-panel.vcf')

mutFilterPON <- function(rawmaf, vcf_add) {
  
  # first filter using the flag in maf column
  if('PON' %in% colnames(rawmaf)) {
    rawmaf <- data.frame(rawmaf[is.na(rawmaf$PON), ])
  }
  rownames(rawmaf) <- 1:nrow(rawmaf)
  
  discard <- c()
  
  # joint useful information for matching
  if (unique(rawmaf$NCBI_Build) == 'GRCh37') {
    ext37 <- rep(NA, nrow(somatic37@fix))
    for (k in 1:nrow(somatic37@fix)) {
      ext37[k] <- paste(somatic37@fix[k, c(1, 2, 4, 5)], collapse = ";")
    }
    
    # combine useful info for matching
    for (i in 1:nrow(rawmaf)) {
      info <- paste(c(str_extract(rawmaf$Chromosome[i], '[:digit:]+'), 
                      rawmaf[i, c(6, 11, 13)]), collapse = ";")
      if (info %in% ext37) {
        discard <- c(discard, i)
      }
    }
  } else if (unique(rawmaf$NCBI_Build) == 'GRCh38') {
    
    ext38 <- rep(NA, nrow(somatic38@fix))
    for (k in 1:nrow(somatic38@fix)) {
      ext38[k] <- paste(somatic38@fix[k, c(1, 2, 4, 5)], collapse = ";")
    }
    
    # combine useful info for matching
    for (i in 1:nrow(rawmaf)) {
      info <- paste(c(str_extract(rawmaf$Chromosome[i], '[:digit:]+'), 
                      rawmaf[i, c(6, 11, 13)]), collapse = ";")
      if (info %in% ext38) {
        discard <- c(discard, i)
      }
    }
    
  } else if (!(unique(rawmaf$NCBI_Build) %in% c('GRCh38', 'GRCh37'))){
    stop('Invaild NCBI build info')
  }else if (length(unique(rawmaf$NCBI_Build)) > 1){
    stop('There are more than one NCBI build version in this dataset')
  }
  
  # discard somatic mutation
  maf <- as.data.frame(rawmaf[-discard, ])
  vcf_final <- as.data.frame(vcf_add[-discard, ])
  if (nrow(maf) == 0) {
    
    message('Note: no mutation left after PON filtering')
    return(list(maf, vcf_final))
  } else {
    rownames(maf) <- 1:nrow(maf)
    rownames(vcf_final) <- 1:nrow(vcf_final)
    return(list(maf, vcf_final))
  }
  
}
