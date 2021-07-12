

vcfsToMAF <- function(vcfs_path, tumorSampleNames, jointMAF = TRUE){
  filenames <- list.files(vcfs_path,
                          pattern="*.vep.vcf", full.names=TRUE)
  if (length(filenames) == 0){
    stop('No VCF file detected under the path you offered.')
  }else{
    dfs <- lapply(filenames, vcfToMAF)
    if (length(dfs) == length(tumorSampleNames)){
      names(dfs) <- tumorSampleNames
    }else{
      warning('Length of objects doesn\'t match')
    }
  }

  if (jointMAF) {
    for (i in 1:length(dfs)){
      colnames(dfs[[i]])[114] <- 'TumorSampleInfo'
    }
  }else{
    return(dfs)
  }

}
