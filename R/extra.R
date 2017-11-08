

readRA <- function(genofile, gform, mum, dad, excSamp=NULL, nClust=3,
                   MAFthres=0.05, MISSthres=0.8, SAMPthres=0.01, BINthres=0, DEPTHthres=6, pvalue=0.05){
  
  if(!is.character(genofile) || length(genofile) != 1)
    stop("File name of RA data set is not a string of length one")
  if(!is.character(gform) || length(gform) != 1 || !(gform %in% c("Tassel","uneak")))
    stop("gform argument must be either 'Tassel' ot 'uneak'")
  
  ### check whether the filtering thresholds have been 
  
  
  ## separate character between reference and alternate allele count
  gsep <- switch(gform, uneak = "|", Tassel = ",")
  
  ## Process the individuals info
  ghead <- scan(genofile, what = "", nlines = 1, sep = "\t")
  ## index the parents
  mumIndx <- which(ghead %in% mum)
  if(length(mumIndx) == 0)
    stop("Mother ID not found in the data set")
  dadIndx <- which(ghead %in% dad)
  if(length(dadIndx) == 0)
    stop("Father ID not found in the data set")
  progIndx <- switch(gform, uneak = ghead[2:(nInd + 1)], Tassel = (1:length(ghead))[-c(1:2,mumIndx,dadIndx)])
  nInd <- length(ghead) - length(c(mumIndx,dadIndx)) - switch(gform, uneak = 6, Tassel = 2)
  indID <- ghead[progIndx]
  ## Read in the data
  # If reference based
  if (gform == "Tassel"){
    genosin <- scan(genofile, skip = 1, sep = "\t", what = c(list(chrom = "", coord = 0), rep(list(""), nInd+length(c(mumIndx,dadIndx)))))
    chrom <- genosin[[1]]
    pos <- genosin[[2]]
    SNP_Names <- paste(genosin[[1]],genosin[[2]],sep="_")
  }
  # Non-reference based 
  if (gform == "uneak"){
    genosin <- scan(genofile, skip = 1, sep = "\t", what = c(list(chrom = ""), rep(list(""), nind), list(hetc1 = 0, hetc2 = 0, acount1 = 0, acount2 = 0, p = 0)))
    SNP_Names <- genosin[[1]]
  }
  ## Compute the number of SNPs
  nSnps <- length(SNP_Names)
  
  ## Compute the read count matrices and genon matrix for the offspring
  depth_Ref <- depth_Alt <- matrix(0, nrow = nInd, ncol = nSnps)
  for (i in 1:nInd){ 
    depths <- strsplit(genosin[[progIndx[i]]], split = gsep, fixed = TRUE)
    depth_Ref[i, ] <- as.numeric(unlist(lapply(depths,function(z) z[1])))
    depth_Alt[i, ] <- as.numeric(unlist(lapply(depths,function(z) z[2])))
  }
  genon <- (depth_Ref > 0) + (depth_Alt == 0)
  genon[which(depth_Ref == 0 & depth_Alt == 0)] <- NA
  if (gform == "uneak") 
    AFrq <- genosin[[length(genosin)]]
  
  ## Check that samples meet the sample treshold
  sampDepth <- rowMeans(depth_Ref + depth_Alt)
  excSamp <- unique(c(excSamp,indID[which(sampDepth < SAMPthres)]))
  
  ## Remove any sample which we don't want
  if(!is.null(excSamp)){
    toRemove <- which(indID %in% excSamp)
    if(length(excSamp) > 0){
      depth_Ref <- depth_Ref[-toRemove,]
      depth_Alt <- depth_Alt[-toRemove,]
      genon <- genon[-toRemove,]
      indID <- indID[-toRemove]
      nInd <- length(indID)
    }
  }
  
  ## Determine the segregation types of the loci
  genon_mum <- depth_mum <- matrix(nrow=length(mumIndx),ncol=nSnps)
  for(i in 1:length(mumIndx)){
    genon_mum[i,] <- unlist(lapply(strsplit(genosin[[mumIndx[i]]],split=","), function(x) switch(1+2*(x[1]==0) + (x[2]==0),1,2,0,NA)))
    depth_mum[i,] <- unlist(lapply(strsplit(genosin[[mumIndx[i]]],split=","), function(x) sum(as.numeric(x))))
  }
  genon_dad <- depth_dad <- matrix(nrow=length(dadIndx),ncol=nSnps)
  for(i in 1:length(dadIndx)){
    genon_dad[i,] <- unlist(lapply(strsplit(genosin[[dadIndx[i]]],split=","), function(x) switch(1+2*(x[1]==0) + (x[2]==0),1,2,0,NA)))
    depth_dad[i,] <- unlist(lapply(strsplit(genosin[[dadIndx[i]]],split=","), function(x) sum(as.numeric(x))))
  }
  
  ## Determine segregation type of each SNP if possible
  config <- unlist(sapply(1:nSnps,function(x){
    x_p = genon_dad[,x]; x_m = genon_mum[,x]
    d_p = depth_dad[,x]; d_m = depth_mum[,x]
    if(sum(d_p)>DEPTHthres & sum(d_p)>DEPTHthres ){
      if(any(x_p==1,na.rm=T) & any(x_m==1,na.rm=T))
        return(1)
      else if(all(x_m==2,na.rm=T) & (any(x_p==1, na.rm=T) | all(x_p %in% c(0,2), na.rm=T)))
        return(2)
      else if(all(x_m==0,na.rm=T) & (any(x_p==1, na.rm=T) | all(x_p %in% c(0,2), na.rm=T)))
        return(3)
      else if(all(x_p==2,na.rm=T) & (any(x_m==1, na.rm=T) | all(x_m %in% c(0,2), na.rm=T)))
        return(4)
      else if(all(x_p==0,na.rm=T) & (any(x_m==1, na.rm=T) | all(x_m %in% c(0,2), na.rm=T)))
        return(5)
      else
        return(NA)
    }
    else return(NA)
  }))
  
  #### Segregation test to determine if the SNPs have been miss-classified
  seg_Dis <- sapply(1:nSnps,function(x){
    if(is.na(config[x]))
      return(NA)
    else{
      d = depth_Ref[,x] + depth_Alt[,x]
      g = genon[,x]
      K = sum(1/2^(d[which(d != 0)])*0.5)/sum(d != 0)
      nAA = sum(g==2, na.rm=T)
      nAB = sum(g==1, na.rm=T)
      nBB = sum(g==0, na.rm=T)
      ## check that there are sufficient data to perform the chisq test
      if(sum(nAA+nAB+nBB)/length(g) <= (1-MISSthres))
        return(NA)
      else if(config[x] == 1){
        exp_prob <- c(0.25 + K,0.5 - 2*K, 0.25 + K)
        ctest <- chisq.test(c(nBB,nAB,nAA), p = exp_prob)
        return(ifelse(ctest$p.value < pvalue, TRUE, FALSE))
      }
      else if(config[x] %in% c(2,4)){
        exp_prob <- c(K, 0.5 - 2*K, 0.5 + K)
        ctest <- chisq.test(c(nBB,nAB,nAA), p = exp_prob)
        return(ifelse(ctest$p.value < pvalue, TRUE, FALSE))
      }
      else if(config[x] %in% c(3,5)){
        exp_prob <- c(0.5 + K, 0.5 - 2*K, K)
        ctest <- chisq.test(c(nBB,nAB,nAA), p = exp_prob)
        return(ifelse(ctest$p.value < pvalue, TRUE, FALSE))
      }
    }
  },simplify = T)
  config[which(seg_Dis)] <- NA
  
  cat("Performing filtering:\n")
  cat("Criteria for retaining SNPs are:\n")
  cat("Minor allele frequency (MAF) >",MAFthres,"\n")
  cat("Percentage of missing genotypes < ",MISSthres*100,"%\n",sep="")
  ## Run the filtering of the progeny SNPs
  MAF <- colMeans(genon, na.rm=T)/2
  MAF <- pmin(MAF,1-MAF)
  miss <- apply(genon,2, function(x) sum(is.na(x))/length(x))
  
  ## Infer geotypes for over SNPs that have passed the MAF and MISS thresholds
  propHeter <- sapply(1:nSnps, function(x) sum(genon[,x] == 1,na.rm=T)/sum(!is.na(genon[,x])))
  toInfer <- (MAF > MAFthres) & (miss < MISSthres) & (propHeter < 0.9) & is.na(config)
  
  seg_Infer <- sapply(1:nSnps, function(x){
    if(!toInfer[x])
      return(NA)
    else{
      d = depth_Ref[,x] + depth_Alt[,x]
      g = genon[,x]
      K = sum(1/2^(d[which(d != 0)])*0.5)/sum(d != 0)
      nAA = sum(g==2, na.rm=T)
      nAB = sum(g==1, na.rm=T)
      nBB = sum(g==0, na.rm=T)
      ## check that there are sufficient data to perform the chisq test
      if(sum(nAA+nAB+nBB)/length(g) <= (1-MISSthres))
        return(NA)
      ## compute chiseq test for both loci types
      exp_prob_BI <- c(0.25 + K,0.5 - 2*K, 0.25 + K)
      exp_prob_SI <- c(K, 0.5 - 2*K, 0.5 + K)
      ctest_BI <- chisq.test(c(nBB,nAB,nAA), p = exp_prob_BI)
      ctest_SI_1 <- chisq.test(c(nBB,nAB,nAA), p = exp_prob_SI)
      ctest_SI_2 <- chisq.test(c(nBB,nAB,nAA), p = rev(exp_prob_SI))
      ## do tests to see if we can infer type
      if( ctest_BI$p.value > pvalue & ctest_SI_1$p.value < pvalue & ctest_SI_2$p.value < pvalue )
        return(1)
      else if ( ctest_BI$p.value < pvalue & ctest_SI_1$p.value > pvalue & ctest_SI_2$p.value < pvalue )
        return(4)
      else if ( ctest_BI$p.value < pvalue & ctest_SI_1$p.value < pvalue & ctest_SI_2$p.value > pvalue )
        return(5)
      else
        return(NA)
    }
  },simplify = T)
  
  indx <- (MAF > MAFthres) & (miss < MISSthres) & ( !is.na(config) | !is.na(seg_Infer) )
  
  ## Determine the segregation of the groups
  group <- list()
  group$BI <- which(config[indx] == 1)
  group$PI <- which(config[indx] %in% c(2,3))
  group$MI <- which(config[indx] %in% c(4,5))
  
  group_infer <- list()
  group_infer$BI <- which(seg_Infer[indx] == 1) 
  group_infer$SI <- which(seg_Infer[indx] %in% c(4,5))
  
  # if(gform == "Tassel" & BINthres > 0){
  #   ## Extract single snp from each tag
  #   oneSNP <- sapply(unique(chrom), function(x){
  #     ind <- which(chrom[indx] == x)
  #     if(length(ind) == 0)
  #       print("next")
  #     else{
  #       g1_diff <- diff(pos.filt[ind])
  #       SNP_bin <- c(0,cumsum(g1_diff > BINthres)) + 1
  #       set.seed(58473+ind[1])
  #       pos <- sapply(unique(SNP_bin), function(y) {
  #         ind2 <- which(SNP_bin == y)
  #         if(length(ind2) > 1)
  #           return(sample(ind2,size=1))
  #         else if(length(ind2) == 1)
  #           return(ind2)
  #       })
  #       return(ind[pos])
  #     }
  #   } )
  # }
  
  # Subset the data
  depth_Ref <- depth_Ref[,indx]
  depth_Alt <- depth_Alt[,indx]
  genon <- genon[,indx]
  chrom <- chrom[indx]
  pos <- pos[indx]
  SNP_Names <- SNP_Names[indx]
  config <- config[indx]
  config_infer <- seg_Infer[indx]
  nSnps <- sum(indx)
  
  cat("\n\nSummary:\n\n")
  cat("Number of SNPs remaining after filtering:",nSnps,"\n")
  cat("Number of progeny:", nInd,"\n")
  cat("Number of SNPs with correct segregation type:", sum(!is.na(config)),"\n")
  cat("Both-informative (BI):", length(group$BI),"\n")
  cat("Maternal-informative (MI):", length(group$MI),"\n")
  cat("Paternal-informative (PI):", length(group$PI),"\n")
  cat("Number of SNPs with inferred segregation type:", sum(!is.na(config_infer)),"\n")
  cat("Both-informative (BI):", length(group_infer$BI),"\n")
  cat("Maternal/Paternal-informative (MI or PI):", length(group_infer$SI),"\n")
  
      
  obj <- list(genon=genon, depth_Ref=depth_Ref,depth_Alt=depth_Alt,chrom=chrom,pos=pos,indID=indID,SNP_Names=SNP_Names,config=config, 
              config_infer=config_infer, group=group, group_infer=group_infer, nInd=nInd, nSnps=nSnps)
  
  ## Add filtering at some stage
  obj$filtering <- NULL
    
  ## Compute the 2-point recombination fraction estimates
  obj <- rf_2pt(obj, nClust=nClust)
  
  ## return an object of the type we want
  return(obj)
}

## needed for foreach loop
comb <- function(...){
  mapply('rbind',...,SIMPLIFY=FALSE)
}

rf_2pt <- function(obj, nClust, inferSNPs = TRUE){
  
  if(length(c(unlist(obj$group), unlist(obj$group_infer))) == 0)
    stop("There are no SNPs in the data set.")
  
  ## Cacluate the number of SNPs
  if(length(obj$group_infer$BI) > 0 & inferSNPs){
    indx_BI <- c(obj$group$BI, obj$group_infer$BI)
    nSnps_BI <- length(indx_BI)
  }
  else{
    indx_BI <- c(obj$group$BI)
    nSnps_BI <- length(indx_BI)
  }
  if(length(obj$group_infer$SI) > 0 & inferSNPs){
    indx_MI <- c(obj$group$MI, obj$group_infer$SI)
    nSnps_MI <- length(indx_MI)
  }
  else{
    indx_MI <- c(obj$group$MI)
    nSnps_MI <- length(indx_MI)
  }
  indx_PI <- c(obj$group$PI)
  nSnps_PI = length(indx_PI)
  
  ## set up the config vector
  config <-  obj$config
  config[which(!is.na(obj$config_infer))] <- obj$config_infer[which(!is.na(obj$config_infer))]
  if(any(is.na(config)))
    stop("There are some missing segregation types in the data.")
  
  ## Set up the Clusters
  cl <- makeCluster(nClust)
  registerDoSNOW(cl)
  
  cat("\nComputing 2-point recombination fraction estimates ...\n")
  cat("Paternal SNPs\n")
  ## Paternal informative SNPs
  rf.PI <- foreach(snp1 = iter(seq(length.out=nSnps_PI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_PI),simplify=F)
    for(snp2 in seq_len(snp1-1)){
      ind = indx_PI[c(snp1,snp2)]
      rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                                    OPGP=list(c(5,5) + 2*(config[ind]==3)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                                    OPGP=list(c(5,6) + 2*(config[ind]==3)), epsilon=NULL)
      rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
      rf[[1]][snp2] <- rf.ind$rf
      rf[[2]][snp2] <- rf.ind$LOD
    }
    return(rf)
  }
  for(i in 1:2){
    rf.PI[[i]][upper.tri(rf.PI[[i]])] <- t(rf.PI[[i]])[upper.tri(rf.PI[[i]])]
  }

  cat("Maternal SNPs\n")
  ## Maternal informative SNPs
  rf.MI <- foreach(snp1 = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_MI),simplify=F)
    for(snp2 in seq_len(snp1-1)){
      ind = indx_MI[c(snp1,snp2)]
      rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                                    OPGP=list(c(9,9) + 2*(config[ind]==5)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                                    OPGP=list(c(9,10) + 2*(config[ind]==5)), epsilon=NULL)
      rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
      rf[[1]][snp2] <- rf.ind$rf
      rf[[2]][snp2] <- rf.ind$LOD
    }
    return(rf)
  }  
  for(i in 1:2){
    rf.MI[[i]][upper.tri(rf.MI[[i]])] <- t(rf.MI[[i]])[upper.tri(rf.MI[[i]])]
  }
  
  cat("Both informative SNPs\n")
  ### Both Informative
  rf.BI <- foreach(snp1 = iter(seq(length.out=nSnps_BI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_BI),simplify=F)
    for(snp2 in seq_len(snp1-1)){
      ind = indx_BI[c(snp1,snp2)]
      temp1 <- rf_est_FS(0.1,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,1)), epsilon=NULL)
      temp2 <- rf_est_FS(0.4,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,1)), epsilon=NULL)
      rf.est1 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
      temp1 <- rf_est_FS(0.1,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,2)), epsilon=NULL)
      temp2 <- rf_est_FS(0.4,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,2)), epsilon=NULL)
      rf.est2 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
      temp1 <- rf_est_FS(0.1,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,4)), epsilon=NULL)
      temp2 <- rf_est_FS(0.4,depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]), OPGP=list(c(1,4)), epsilon=NULL)
      rf.est4 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
      rf.ind <- switch(which.min(c(rf.est1$loglik, rf.est2$loglik, rf.est4$loglik) ),
                         rf.est1,rf.est2,rf.est4)
      rf[[1]][snp2] <- rf.ind$rf
      rf[[2]][snp2] <- rf.ind$LOD
    }
    return(rf)
  }
  for(i in 1:2){
    rf.BI[[i]][upper.tri(rf.BI[[i]])] <- t(rf.BI[[i]])[upper.tri(rf.BI[[i]])]
  }
  
  cat("Paternal information vs Both informative\n")
  ## Paternal and Informative SNPs
  rf.PI.BI <- foreach(snp.ps = iter(seq(length.out=nSnps_PI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_BI),simplify=F)
    for(snp.bi in 1:nSnps_BI){
      ind <- c(indx_PI[snp.ps],indx_BI[snp.bi])
      rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                           OPGP=list(c(5,1) + 2*c(config[ind[1]]==3,0)), epsilon=NULL )
      rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                                    OPGP=list(c(5,2) + 2*c(config[ind[1]]==3,0)), epsilon=NULL )
      rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
      rf[[1]][snp.bi] <- rf.ind$rf
      rf[[2]][snp.bi] <- rf.ind$LOD
    }
    return(rf)
  }
  
  cat("Maternal information vs Both informative\n")
  ## Maternal and Informative SNPs
  rf.MI.BI <- foreach(snp.ms = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_BI),simplify=F)
    for(snp.bi in 1:nSnps_BI){
      ind <- c(indx_MI[snp.ms],indx_BI[snp.bi])
      rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                           OPGP=list(c(9,1) + 2*c(config[ind[1]]==5,0)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                           OPGP=list(c(9,3) + 2*c(config[ind[1]]==5,0)), epsilon=NULL)
      rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
      rf[[1]][snp.bi] <- rf.ind$rf
      rf[[2]][snp.bi] <- rf.ind$LOD
    }
    return(rf)
  }
  
  ## For the non-informative computations
  ## Really done so that we can check that there is no miss identification of the group
  cat("Maternal information vs Paternal informative\n")
  rf.MI.PI <- foreach(snp.ms = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_PI),simplify=F)
    for(snp.pi in 1:nSnps_PI){
      ind <- c(indx_MI[snp.ms],indx_PI[snp.pi])
      rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                           OPGP=list(c(9,9) + 2*(config[ind] %in% c(3,5))), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                           OPGP=list(c(9,10) + 2*(config[ind] %in% c(3,5))), epsilon=NULL)
      rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
      rf[[1]][snp.pi] <- rf.ind$rf
      rf[[2]][snp.pi] <- rf.ind$LOD
    }
    return(rf)
  }
  stopCluster(cl) 
  
  ## Build the rf and LOD matrices
  origOrder <- order(c(indx_BI,indx_PI,indx_MI))
  rf.mat <- rbind(cbind(rf.BI[[1]],t(rf.PI.BI[[1]]),t(rf.MI.BI[[1]])),
                        cbind(rf.PI.BI[[1]], rf.PI[[1]], t(rf.MI.PI[[1]])),
                        cbind(rf.MI.BI[[1]], rf.MI.PI[[1]], rf.MI[[1]]))[origOrder,origOrder]
  LOD.mat <- rbind(cbind(rf.BI[[2]],t(rf.PI.BI[[2]]),t(rf.MI.BI[[2]])),
                        cbind(rf.PI.BI[[2]], rf.PI[[2]], t(rf.MI.PI[[2]])),
                        cbind(rf.MI.BI[[2]], rf.MI.PI[[2]], rf.MI[[2]]))[origOrder,origOrder]
  obj$rf <- rf.mat
  obj$LOD <- LOD.mat
  return(obj)
}


### Function for applying additional filters to SNPs (only works after running through readRA)
excludeSNPs <- function(obj, criteria, name){
  if(!is.null(names(obj$filter)) && any(names(obj$filter) == name))
    stop("Filter criteria already exists")
  else
    obj$filter[[name]] <- !criteria(obj)
  return(obj)
}

### Function for creating Linkage groups
createGroups <- function(obj, parent, LOD=10, nComp=10){
  
  ## Initialize the LGs
  LGs <- list()
  
  ## check that the parent argument is correct
  if(!is.character(parent) || length(parent) != 1 || 
     !(parent %in% c("maternal","paternal")))
    stop("parent argument is not a string of of length one or is incorrect:
         Please select one of the following:
         maternal: Only MI SNPs
         paternal: Only PI SNPs")
  
  if(parent == "maternal")
    unmapped <- sort(c(obj$group$MI))
  else if(parent == "paternal")
    unmapped <- sort(c(obj$group$PI))
  
  ## Run algorithm for generating the linkage groups
  finish = FALSE
  while(!finish){
    newLG <- unmapped[sort(which(obj$LOD[unmapped,unmapped]==max(obj$LOD[unmapped,unmapped],na.rm=T),arr.ind=T)[1,])]
    unmapped <- unmapped[-which(unmapped%in%newLG)]
    compLOD <- obj$LOD[newLG,unmapped]
    meanLOD <- apply(compLOD,2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
    maxMeanLOD <- max(meanLOD)
    while(maxMeanLOD > LOD ){
      newSNP <- which(meanLOD == maxMeanLOD)
      newLG <- c(newLG,unmapped[newSNP])
      unmapped <- unmapped[-newSNP]
      compLOD <- obj$LOD[newLG,unmapped]
      meanLOD <- apply(compLOD,2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
      maxMeanLOD <- max(meanLOD)
    }
    meanLOD <- apply(obj$LOD[unmapped,unmapped],2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
    finish <- !any(meanLOD > LOD)
    LGs <- c(LGs,list(newLG))
  }
  
  #obj$LGs <- LGs
  return(LGs)
}

addSNPs <- function(obj, LG.list, LOD=10, nComp=10){
  
  if(!is.list(LG.list))
    stop("The Linkage group object needs to be a list")
  if(length(LG.list) == 0)
    stop("There are no linkage groups. Please use 'createGroups' function to create the groups")
  nLGs <- length(LG.list)
  
  ## Find the unmapped loci
  unmapped <- sort(unlist(obj$group_infer$SI))
  if(length(unmapped) == 0)
    stop("There are no SNPs remaining that are unmapped")
  
  ## Run algorithm for generating the linkage groups
  noneMapped = FALSE
  while(!noneMapped){
    noneMapped = TRUE
    ## check that there are still SNPs remaining that need to be mapped
    if(length(unmapped) == 0)
      return(LG.list)
    ## run the algorithm to map the SNPs
    else{
      for(snp in unmapped){
        LODvalue = numeric(nLGs)
        for(lg in 1:nLGs)
          LODvalue[lg] <- mean(sort(obj$LOD[snp,LG.list[[lg]]],decreasing=T)[1:nComp],na.rm=T)
        if(max(LODvalue) >= LOD & sort(LODvalue, decreasing = T)[2] < LOD){
          newLG <- which.max(LODvalue)
          LG.list[[newLG]] <- c(LG.list[[newLG]], snp)
          unmapped <- unmapped[-which(unmapped == snp)]
          noneMapped = FALSE
        }
      }
    }
  }
  return(LG.list)
}

addBIsnps <- function(obj, LG.list, LOD=10, nComp=10){
  
  if(!is.list(LG.list))
    stop("The Linkage group object needs to be a list")
  if(length(LG.list) == 0)
    stop("There are no linkage groups. Please use 'createGroups' function to create the groups")
  nLGs <- length(LG.list)
  
  LG.list.new <- LG.list
  
  ## Find the unmapped loci
  unmapped <- sort(unlist(obj$group$BI,obj$group_infer$BI))
  if(length(unmapped) == 0)
    stop("There are no SNPs remaining that are unmapped")
  
  ## Run algorithm for generating the linkage groups
  noneMapped = FALSE
  while(!noneMapped){
    noneMapped = TRUE
    ## check that there are still SNPs remaining that need to be mapped
    if(length(unmapped) == 0)
      return(LG.list)
    ## run the algorithm to map the SNPs
    else{
      for(snp in unmapped){
        LODvalue = numeric(nLGs)
        for(lg in 1:nLGs)
          LODvalue[lg] <- mean(sort(obj$LOD[snp,LG.list[[lg]]],decreasing=T)[1:nComp],na.rm=T)
        if(max(LODvalue) >= LOD){
          newLG <- which.max(LODvalue)
          LG.list.new[[newLG]] <- c(LG.list.new[[newLG]], snp)
          unmapped <- unmapped[-which(unmapped == snp)]
          noneMapped = FALSE
        }
      }
    }
  }
  return(LG.list.new)
}

    
#### Function for plotting linkage groups (or a single linkage group)
plotLGs <- function(obj, LG.list, filename=NULL, names=NULL, chrS=2, lmai=2, chrom=T){
  
  rf.mat <- obj$rf
  b <- ncol(rf.mat) + 1
  if(chrom)
    chrom.ind <- unlist(lapply(LG.list, function(x) c(x,b)))[-length(unlist(LG.list))+length(LG.list)]
  else
    chrom.ind <- 1:ncol(rf.mat)
  
  ## Subset the matrix 
  rf.mat <- cbind(rf.mat,rep(0,b-1))
  rf.mat <- rbind(rf.mat,rep(0,b))       
  rf.mat <- rf.mat[chrom.ind,chrom.ind]
  ## work out where the breaks are
  breaks <- which(chrom.ind==b)
  npixels <- length(chrom.ind)
  if(chrom){
    if(!is.null(filename))
      png(filename,width=npixels+72*lmai,height=npixels,res=72)
    par(xaxt='n',yaxt='n',mai=c(0,lmai,0,0),bty='n',ann=F)
    image(1:npixels,1:npixels,rf.mat,zlim=c(0,0.5),col=heat.colors(100))
    if(is.null(names))
      mtext(paste("LG",1:length(LG.list),"  "), 
                  at=floor(apply(cbind(c(0,breaks),c(breaks,npixels)),1,median)),side=2, line=0,cex=chrS,las=1)
    else
      mtext(names, at=floor(apply(cbind(c(0,breaks),c(breaks,npixels)),1,median)),side=2, line=0,cex=chrS,las=1)
    abline(h=breaks)
    abline(v=breaks)
    if(!is.null(filename))
      dev.off()
  }
  else{
    npixels <- length(chrom.ind)
    if(!is.null(filename))
      png(filename,width=npixels,height=npixels)
    par(xaxt='n',yaxt='n',mar=c(0,0,0,0),bty='n',ann=F)
    image(1:npixels,1:npixels,rf.mat,zlim=c(0,0.5),col=heat.colors(100))
    abline(h=breaks)
    abline(v=breaks)
    if(!is.null(filename))
      dev.off()
  }
}


mergeGroups <- function(obj, LOD=3){
  print("To be implemented")
}


#### Function for doing the Ordering
orderLG <- function(obj, LG.list, sigma=10){
  if(!is.list(LG.list))
    stop("Linkage group object is not a list")
  ## Define new list for ordered SNPs
  LG.mat.ord <- list()
  # iterate over the linkage groups
  for(lg in 1:length(LG.list)){
    snpInd <- LG.list[[lg]]
    ## optimization based on initial ordering
    D1 <- as.dist(obj$rf[snpInd,snpInd])
    out1 <- seriate.distLD(D1,method="SPIN_NH_LD",control=list(sigma=sigma,verbose=T))
    ## Try two other random orders
    set.seed(58392+lg*9/3)
    order2 <- sample(1:length(snpInd))
    if(length(snpInd) > 20){
      nSect <- floor(length(snpInd)/10)
      order3 <- order(as.numeric(as.character(cut(1:length(snpInd),breaks=nSect,labels=sample(1:nSect)))))
    }
    else
      order3 <- sample(1:length(snpInd))
    D2 <- as.dist(obj$rf[snpInd[order2],snpInd[order2]])
    D3 <- as.dist(obj$rf[snpInd[order3],snpInd[order3]])
    out2 <- seriate.distLD(D2,method="SPIN_NH_LD",control=list(sigma=sigma,verbose=T))
    out3 <- seriate.distLD(D3,method="SPIN_NH_LD",control=list(sigma=sigma,verbose=T))
    ## Find out which has the best ordering and output that order
    bestOut <- switch(which.min(c(out1[[2]],out2[[2]],out3[[2]])),1,2,3)
    if(bestOut == 1) {
      LG.mat.ord[[lg]] <- snpInd[get_order(out1[[1]])]
    } else if(bestOut == 2){
      LG.mat.ord[[lg]] <- snpInd[order2][get_order(out2[[1]])]
    } else if(bestOut == 3){
      LG.mat.ord[[lg]] <- snpInd[order3][get_order(out3[[1]])]
    }
  }
  return(LG.mat.ord)
}

