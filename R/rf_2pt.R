
## function needed for foreach loop
comb <- function(...){
  mapply('rbind',...,SIMPLIFY=FALSE)
}

### Function for computing the pairwise RF in a single full-sib family.
rf_2pt_single <- function(depth_Ref, depth_Alt, config, config_infer, group, group_infer, inferSNPs, nClust){
  
  if(length(c(unlist(group), unlist(group_infer))) == 0)
    stop("There are no SNPs in the data set.")
  
  ## Calcuate the number of SNPs
  if(length(group_infer$BI) > 0 & inferSNPs){
    indx_BI <- c(group$BI, group_infer$BI)
    nSnps_BI <- length(indx_BI)
  }
  else{
    indx_BI <- c(group$BI)
    nSnps_BI <- length(indx_BI)
  }
  if(length(group_infer$SI) > 0 & inferSNPs){
    indx_MI <- c(group$MI, group_infer$SI)
    nSnps_MI <- length(indx_MI)
  }
  else{
    indx_MI <- c(group$MI)
    nSnps_MI <- length(indx_MI)
  }
  indx_PI <- c(group$PI)
  nSnps_PI = length(indx_PI)
  
  ## set up the config vector
  if(inferSNPs)
    config[which(!is.na(config_infer))] <- config_infer[which(!is.na(config_infer))]
  ## Check that there are no missing configs in the data
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
      rf.est1 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
                           OPGP=list(c(5,5) + 2*(config[ind]==3)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
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
      rf.est1 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
                           OPGP=list(c(9,9) + 2*(config[ind]==5)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
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
      temp1 <- rf_est_FS(0.1,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,1)), epsilon=NULL)
      temp2 <- rf_est_FS(0.4,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,1)), epsilon=NULL)
      rf.est1 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
      temp1 <- rf_est_FS(0.1,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,2)), epsilon=NULL)
      temp2 <- rf_est_FS(0.4,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,2)), epsilon=NULL)
      rf.est2 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
      temp1 <- rf_est_FS(0.1,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,4)), epsilon=NULL)
      temp2 <- rf_est_FS(0.4,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,4)), epsilon=NULL)
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
      rf.est1 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
                           OPGP=list(c(5,1) + 2*c(config[ind[1]]==3,0)), epsilon=NULL )
      rf.est2 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
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
      rf.est1 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
                           OPGP=list(c(9,1) + 2*c(config[ind[1]]==5,0)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
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
      rf.est1 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
                           OPGP=list(c(9,9) + 2*(config[ind] %in% c(3,5))), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
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
  return(list(rf = rf.mat,LOD = LOD.mat))
}



### Function for computing the pairwise RF in a multiple-sib family.
rf_2pt_multi <- function(depth_Ref, depth_Alt, config, config_infer = NULL, group, group_infer, inferSNPs){
  
  if(length(c(unlist(group), unlist(group_infer))) == 0)
    stop("There are no SNPs in the data set.")
  
  ## Calcuate the number of SNPs
  if(length(group_infer$BI) > 0 & inferSNPs){
    indx_BI <- c(group$BI, group_infer$BI)
    nSnps_BI <- length(indx_BI)
  }
  else{
    indx_BI <- c(group$BI)
    nSnps_BI <- length(indx_BI)
  }
  if(length(group_infer$SI) > 0 & inferSNPs){
    indx_MI <- c(group$MI, group_infer$SI)
    nSnps_MI <- length(indx_MI)
  }
  else{
    indx_MI <- c(group$MI)
    nSnps_MI <- length(indx_MI)
  }
  indx_PI <- c(group$PI)
  nSnps_PI = length(indx_PI)
  
  ## set up the config vector
  config <-  config
  if(inferSNPs)
    config[which(!is.na(config_infer))] <- config_infer[which(!is.na(config_infer))]
  ## Check that there are no missing configs in the data
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
      rf.est1 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
                           OPGP=list(c(5,5) + 2*(config[ind]==3)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
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
      rf.est1 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
                           OPGP=list(c(9,9) + 2*(config[ind]==5)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
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
      temp1 <- rf_est_FS(0.1,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,1)), epsilon=NULL)
      temp2 <- rf_est_FS(0.4,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,1)), epsilon=NULL)
      rf.est1 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
      temp1 <- rf_est_FS(0.1,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,2)), epsilon=NULL)
      temp2 <- rf_est_FS(0.4,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,2)), epsilon=NULL)
      rf.est2 <- switch(which.min(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
      temp1 <- rf_est_FS(0.1,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,4)), epsilon=NULL)
      temp2 <- rf_est_FS(0.4,depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]), OPGP=list(c(1,4)), epsilon=NULL)
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
      rf.est1 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
                           OPGP=list(c(5,1) + 2*c(config[ind[1]]==3,0)), epsilon=NULL )
      rf.est2 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
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
      rf.est1 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
                           OPGP=list(c(9,1) + 2*c(config[ind[1]]==5,0)), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
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
      rf.est1 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
                           OPGP=list(c(9,9) + 2*(config[ind] %in% c(3,5))), epsilon=NULL)
      rf.est2 <- rf_est_FS(depth_Ref=list(depth_Ref[,ind]),depth_Alt=list(depth_Alt[,ind]),
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
  return(list(rf.mat,LOD.mat))
}


