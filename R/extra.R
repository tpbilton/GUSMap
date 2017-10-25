

readRA <- function(genofile, gform, mum, dad, excSamp=NULL, nClust=3,
                   MAFthres=0.01, MISSthres=0.5, SAMPthres=0.01, BINthres=0, DEPTHthres=0, pvalue=0.2){
  
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
      else if(all(x_m==2,na.rm=T) & any(x_p==1, na.rm=T))
        return(2)
      else if(all(x_m==0,na.rm=T) & any(x_p==1, na.rm=T))
        return(3)
      else if(all(x_p==2,na.rm=T) & any(x_m==1,na.rm=T))
        return(4)
      else if(all(x_p==0,na.rm=T) & any(x_m==1,na.rm=T))
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
  indx <- (MAF > MAFthres) & (miss < MISSthres) & !is.na(config)
  
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
  nSnps <- sum(indx)

  cat("Number of SNPs remaining after filtering:",nSnps,"\n")
  cat("Number of progeny:", nInd,"\n")
  
  ## Determine the segregation of the groups
  group <- list()
  group$BI <- which(config == 1)
  group$PI <- which(config %in% c(2,3))
  group$MI <- which(config %in% c(4,5))
  
  obj <- list(genon=genon, depth_Ref=depth_Ref,depth_Alt=depth_Alt,chrom=chrom,pos=pos,indID=indID,SNP_Names=SNP_Names,config=config,
              group=group, nInd=nInd, nSnps=nSnps)
  
  ## Add filtering at some stage
  obj$filtering <- NULL
    
  ## Compute the 2-point recombination fraction estimates
  obj <- rf_2pt(obj, nClust=nClust)
  
  ## return an object of the type we want
  return(obj)
}

#dataG <- readRA(genofile,"Tassel",mum,dad)

## needed for foreach loop
comb <- function(...){
  mapply('rbind',...,SIMPLIFY=FALSE)
}

rf_2pt <- function(obj, nClust){
  
  ## Cacluate the number of SNPs
  nSnps_BI = length(obj$group$BI)
  nSnps_PI = length(obj$group$PI)
  nSnps_MI = length(obj$group$MI)
  
  ## Set up the Clusters
  cl <- makeCluster(nClust)
  registerDoSNOW(cl)
  
  print("Paternal SNPs\n")
  ## Paternal informative SNPs
  rf.PI <- foreach(snp1 = iter(seq(length.out=nSnps_PI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_PI),simplify=F)
    for(snp2 in seq_len(snp1-1)){
      ind = obj$group$PI[c(snp1,snp2)]
      rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                                    OPGP=list(c(5,5) + 2*(obj$config[ind]==3)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                                    OPGP=list(c(5,6) + 2*(obj$config[ind]==3)), epsilon=NULL)
      rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
      rf[[1]][snp2] <- rf.ind$rf
      rf[[2]][snp2] <- rf.ind$LOD
    }
    return(rf)
  }
  for(i in 1:2){
    rf.PI[[i]][upper.tri(rf.PI[[i]])] <- t(rf.PI[[i]])[upper.tri(rf.PI[[i]])]
  }

  print("Maternal SNPs\n")
  ## Maternal informative SNPs
  rf.MI <- foreach(snp1 = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_MI),simplify=F)
    for(snp2 in seq_len(snp1-1)){
      ind = obj$group$MI[c(snp1,snp2)]
      rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                                    OPGP=list(c(9,9) + 2*(obj$config[ind]==5)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                                    OPGP=list(c(9,10) + 2*(obj$config[ind]==5)), epsilon=NULL)
      rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
      rf[[1]][snp2] <- rf.ind$rf
      rf[[2]][snp2] <- rf.ind$LOD
    }
    return(rf)
  }  
  for(i in 1:2){
    rf.MI[[i]][upper.tri(rf.MI[[i]])] <- t(rf.MI[[i]])[upper.tri(rf.MI[[i]])]
  }
  
  print("Both informative SNPs\n")
  ### Both Informative
  rf.BI <- foreach(snp1 = iter(seq(length.out=nSnps_BI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_BI),simplify=F)
    for(snp2 in seq_len(snp1-1)){
      ind = obj$group$BI[c(snp1,snp2)]
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
  
  print("Paternal information vs Both informative\n")
  ## Paternal and Informative SNPs
  rf.PI.BI <- foreach(snp.ps = iter(seq(length.out=nSnps_PI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_BI),simplify=F)
    for(snp.bi in 1:nSnps_BI){
      ind <- c(obj$group$PI[snp.ps],obj$group$BI[snp.bi])
      rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                           OPGP=list(c(5,1) + 2*c(obj$config[ind[1]]==3,0)), epsilon=NULL )
      rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                                    OPGP=list(c(5,2) + 2*c(obj$config[ind[1]]==3,0)), epsilon=NULL )
      rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
      rf[[1]][snp.bi] <- rf.ind$rf
      rf[[2]][snp.bi] <- rf.ind$LOD
    }
    return(rf)
  }
  
  print("Maternal information vs Both informative\n")
  ## Maternal and Informative SNPs
  rf.MI.BI <- foreach(snp.ms = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_BI),simplify=F)
    for(snp.bi in 1:nSnps_BI){
      ind <- c(obj$group$MI[snp.ms],obj$group$BI[snp.bi])
      rf.est1 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                           OPGP=list(c(9,1) + 2*c(obj$config[ind[1]]==5,0)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(obj$depth_Ref[,ind]),depth_Alt=list(obj$depth_Alt[,ind]),
                           OPGP=list(c(9,3) + 2*c(obj$config[ind[1]]==5,0)), epsilon=NULL)
      rf.ind <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), rf.est1, rf.est2)
      rf[[1]][snp.bi] <- rf.ind$rf
      rf[[2]][snp.bi] <- rf.ind$LOD
    }
    return(rf)
  }
  stopCluster(cl) 
  
  ## For the non-informative computations
  rf.MI.PI <- replicate(2,matrix(NA,nrow=nSnps_MI,ncol=nSnps_PI), simplify=F)
  
  ## Build the rf and LOD matrices
  origOrder <- order(c(obj$group$BI,obj$group$PI,obj$group$MI))
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



createGroups <- function(obj, parent, LOD=10){
  
  ## Initialize the LGs
  LGs <- list()
  
  ## Start with the maternal side
  if(parent=="maternal")
    unmapped <- sort(c(obj$group$MI))
  else(parent=="paternal")
    unmapped <- sort(c(obj$group$PI))
  
  ## Run the process of computing the linkage groups
  finish = FALSE
  while(!finish){
    newLG <- unmapped[sort(which(obj$LOD[unmapped,unmapped]==max(obj$LOD[unmapped,unmapped],na.rm=T),arr.ind=T)[1,])]
    unmapped <- unmapped[-which(unmapped%in%newLG)]
    compLOD <- obj$LOD[newLG,unmapped]
    maxLOD <- max(compLOD)
    while(maxLOD > LOD ){
      newSNP <- matrix(which(compLOD == maxLOD,arr.ind=T),ncol=2)[1,]
      if(newSNP[1]<length(newLG))
        newLG <- c(newLG[1:newSNP[1]],unmapped[newSNP[2]],newLG[(newSNP[1]+1):length(newLG)])
      else
        newLG <- c(newLG,unmapped[newSNP[2]])
      unmapped <- unmapped[-newSNP[2]]
      compLOD <- obj$LOD[newLG,unmapped]
      maxLOD <- max(compLOD)
    }
    finish <- !any(obj$LOD[unmapped,unmapped] > LOD)
    LGs <- c(LGs,list(newLG))
  }
  
  #obj$LGs <- LGs
  return(LGs)
}


mergeGroups <- function(obj, LOD=3){
  print("To be implemented")
}




