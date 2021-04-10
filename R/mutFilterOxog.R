library(pracma)


mutFilterOxog <- function(maf, outFile, globalPoxoG = 0.96, lod0Thresh = -1,
                          artifactThreRate = 0.01, oxoQP1 = 1e8, oxoQP2 = 0){
  
  # check for required columns in maf
  reqCols <- c('Chromosome', 'Start_Position', 'End_Position', 'Variant_Type', 
               'Reference_Allele', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode', 
               'Matched_Norm_Sample_Barcode', 'i_t_ALT_F1R2', 'i_t_ALT_F2R1', 
               'i_t_REF_F1R2', 'i_t_REF_F2R1', 'i_t_Foxog', 'i_picard_oxoQ')
  if (!(all(reqCols %in% colnames(maf)))){
    stop('One or more than one required column is missing in the maf data.')
  }
  
  ## generate data randomly
  # maf <- maf[which(maf$Variant_Type == 'SNP'), ]
  # maf <- cbind(maf, i_t_ALT_F1R2 = round(runif(98, 0, min(maf$t_alt_count))), 
               # i_t_ALT_F2R1 = round(runif(98, 0, 20)), 
               # i_t_REF_F1R2 = round(runif(98, 0, 20)), 
               # i_t_REF_F2R1 = round(runif(98, 0, 20)), 
               # i_t_Foxog = 0, ref_context = 'ATCGGA')
  # for (i in 1:nrow(maf)) {
    # if (maf$Reference_Allele[i] == 'C') {
      # maf$i_t_Foxog[i] <- maf$i_t_ALT_F2R1[i]/(maf$i_t_ALT_F1R2[i] + 
                                                 # maf$i_t_ALT_F2R1[i])
    # }else if(maf$Reference_Allele[i] == 'A'){
      # maf$i_t_Foxog[i] <- maf$i_t_ALT_F2R1[i]/(maf$i_t_ALT_F1R2[i] + 
                                                 # maf$i_t_ALT_F2R1[i])
    # }else if(maf$Reference_Allele[i] == 'G'){
      # maf$i_t_Foxog[i] <- maf$i_t_ALT_F1R2[i]/(maf$i_t_ALT_F1R2[i] + 
                                                 # maf$i_t_ALT_F2R1[i])
    # }else{
      # maf$i_t_Foxog[i] <- maf$i_t_ALT_F1R2[i]/(maf$i_t_ALT_F1R2[i] + 
                                                 # maf$i_t_ALT_F2R1[i])
    # }
  # }
  
  maf <- cbind(maf, i_picard_oxoQ = round(maf$i_t_Foxog, 2))
  
  # Target rate for artifacts that escape filtering
  fdrThresh <- artifactThresholdRate;  
  
  # Measured p used in the B_A(ac,p) distribution in the binomial mixture model
  PoxoG <- globalPoxoG
  
  # Alt Allele Counts to use
  acs <- 3:50
  
  # Null distribution (real mutation) binomial parameter p
  p0 <- 0.5
  
  # pre-processing before filtration
  ## check for tumor-normal pair
  if(nrow(unique(maf[, c('Tumor_Sample_Barcode', 
                         'Matched_Norm_Sample_Barcode')])) != 1){
    stop('More than one tumor-normal sample pair was detected, please check your
         input maf file.')
  }
  
  ## extract and sort maf columns
  mafOxo <- as.data.frame(matrix(ncol = 7, nrow = nrow(maf)))
  colnames(mafOxo) <- c('foxog', 'alt_read_count', 'ref_read_count', 
                        'f1r2', 'f2r1', 'isArtifactMode', 'picard_oxoQ')
  mafOxo$foxog <- maf$i_t_Foxog
  mafOxo$alt_read_count <- maf$i_t_ALT_F1R2 + maf$i_t_ALT_F2R1
  mafOxo$ref_read_count <- maf$i_t_REF_F1R2 + maf$i_t_REF_F2R1
  mafOxo$f1r2 <- maf$i_t_ALT_F1R2
  mafOxo$f2r1 <- maf$i_t_ALT_F2R1
  mafOxo$picard_oxoQ <- maf$i_picard_oxoQ
  for (i in 1:nrow(mafOxo)){
    if ((maf$Reference_Allele[i] == 'C') &(maf$Tumor_Seq_Allele2[i] == 'A') |
        (maf$Reference_Allele[i] == 'G') &(maf$Tumor_Seq_Allele2[i] == 'T')){
      mafOxo$isArtifactMode <- 1
    }else{
      mafOxo$isArtifactMode <- 0
    }
  }
  
  ## recalculate foxog if it equals to -1
  if (any(mafOxo$foxog == -1)){
    mafOxo <- cbind(mafOxo, CA = (maf$Reference_Allele %in% c('C', 'A')))
    CAindex <- rownames(mafOxo[which((mafOxo$foxog == -1) & (mafOxo$CA)), ])
    mafOxo[CAindex, 'foxog'] <- mafOxo[CAindex, 'f2r1'] / 
      mafOxo[CAindex, 'alt_read_count']
    GTindex <- rownames(mafOxo[which((mafOxo$foxog == -1) & !(mafOxo$CA)), ])
    mafOxo[GTindex, 'foxog'] <- mafOxo[GTindex, 'f1r2'] / 
      mafOxo[GTindex, 'alt_read_count']
    mafOxo <- mafOxo[, 1:7]
  }
  
  ## fQ calculation
  oxoQ <- mafOxo$picard_oxoQ[1]
  fQ <- 1/(1+exp(oxoQP2*(oxoQ-oxoQP1)))
  if (oxoQP2 == 0){
    fQ <- 1
  }
  
  # Create a cut using False Discovery Rate 
  ## build FDR frame
  fdrFrame <- as.data.frame(matrix(ncol = 6, nrow = nrow(mafOxo)))
  colnames(fdrFrame) <- c('NALT', 'iART', 'NART', 'pox', 'qox', 'cut')
  fdrFrame$NALT <- mafOxo$alt_read_count
  fdrFrame$iART <- mafOxo$isArtifactMode
  fdrFrame$NART <- round(mafOxo$foxog * mafOxo$alt_read_count)
  fdrFrame$pox <- pbinom(fdrFrame$NART, fdrFrame$NALT, PoxoG)
  fdrFrame[which(fdrFrame$iART == 0), 'pox'] <- 0
  fdrFrame$qox <- calFDR(fdrFrame$pox)
  fdrFrame$cut <- (fdrFrame$qox > fdrThresh)
  
  ## correct for noise floor based on oxoQ estimate
  if (fQ < 1){
    Ncut <- round(sum(fdrFrame$cut)*fQ)
    q <- sort(fdrFrame$pox, decreasing = TRUE)
    i <- order(fdrFrame$pox, decreasing = TRUE)
    #NTOT=length(X.pox);
    #rank(i,1)=(1:NTOT)
    # refine cut selection
    fdrFrame$cut <- fdrFrame$cut & (i <= Ncut)
    
    if (sum(fdrFrame$cut) > 0) {
      # min pox of cut event
      poxmin <- min(fdrFrame$pox(fdrFrame$cut))
      # any pox >= poxmin is cut 
      fdrFrame$cut <- (fdrFrame$cut >= poxmin)
    }
  }
  
  # [Noxo,Nmut,alf,alfci,NoxoCI, lod0] = oxogBinomialMixture(M_A.alt_read_count, round(M_A.foxog .* M_A.alt_read_count), acs, PoxoG, p0);
  
  # Filter
  
  
}


calFDR <- function(p){
  sp <- sort(p)
  ord <- order(p)
  fdr <- (sp * length(p))/1:length(p)
  fdr[fdr > 1] <- 1
  
  fdr <- c(fdr, 1)
  for (i in 1:(length(fdr) - 1)){
    fdr[i] <- min(fdr[i], fdr[i+1])
  }
  fdr <- fdr[-length(fdr)]
  fdr[ord] <- fdr
  return(fdr)
}


oxogBinomialMixture <- function(NALT, NART, AC=1:100, 
                                PoxoG=0.96, p0=0.5, alphaci=0.05){
  
  # drop in negative log likelihood for alphaci (natural log_e)
  dLL <- qnorm(1 - alphaci/2, 0, 1)^2/2
  alpha <- seq(0,1,0.001) 
  #Na=length(alpha);
  #NA=length(AC);
  
  alf <- matrix(nrow = length(AC), ncol = 1)
  alfci <- matrix(nrow = length(AC), ncol = 2)
  Nm <- alf
  cilim <- log10(exp(1))/2
  Nox1 <- alf
  Nox1ci <- alfci
  vlow <- alf
  vhigh <- alf
  
  # initialize lod0
  lod0 <- 0
  
  # define NLL and NLL0 functions
  NLL <- function(a){
    return(sum(-log10(ppois(n, N1*(1-a)*n1+N1*a*n2))))
  }
  
  for (i in 1:length(AC)){
    k <- NALT == i
    nox <- 0:i
    
    if (sum(k) >= 1){
      NALT1 <- NALT[k]
      NART1 <- NART[k]
      N1 <- length(k)
      n1 <- dbinom(nox, i, p0)
      n2 <- dbinom(nox, i, PoxoG)
    }
    
    n <- hist.default(NART1, breaks = c(-Inf, seq(0.5, i-0.5, 1), Inf), 
                      plot = FALSE)$counts
    alf[i] <- fminbnd(NLL, 0, 1)$xmin
    fval0 <- fminbnd(NLL, 0, 1)$fmin
    
    # lod0(i) is log like of fitted alpha - log like of alpha==0
    lod0 <- lod0 - NLL(alf[i]) + NLL(0)
    
    # reduce search region 
    av <- c(0, 10^(-8:-3), seq(0.1,1,0.1))
    fv <- rep(0, length(av))
    for (j in 1:length(av)){
      fv[j] <- NLL(av[j]) - (fval0+dLL)
    }
    
    av <- av[abs(fv) < 1000]
    fv <- fv[abs(fv) < 1000]
  }
}
