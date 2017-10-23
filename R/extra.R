

readRA <- function(genofile, gform, mum, dad, nClust=3){
  
  ## separate character between reference and alternate allele count
  gsep <- switch(gform, uneak = "|", Tassel = ",")
  
  ## Process the individuals info
  ghead <- scan(genofile, what = "", nlines = 1, sep = "\t")
  ## index the parents
  mumIndx <- which(ghead %in% mum)
  dadIndx <- which(ghead %in% dad)
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
  
  ## Determine the segregation types of the loci
  genon_mum <- matrix(nrow=length(mumIndx),ncol=nSnps)
  for(i in 1:length(mumIndx)){
    genon_mum[i,] <- unlist(lapply(strsplit(genosin[[mumIndx[i]]],split=","), function(x) switch(1+2*(x[1]==0) + (x[2]==0),1,2,0,NA)))
  }
  genon_dad <- matrix(nrow=length(dadIndx),ncol=nSnps)
  for(i in 1:length(dadIndx)){
    genon_dad[i,] <- unlist(lapply(strsplit(genosin[[dadIndx[i]]],split=","), function(x) switch(1+2*(x[1]==0) + (x[2]==0),1,2,0,NA)))
  }
  
  config <- unlist(sapply(1:nSnps,function(x){
    x_p = genon_dad[,x]; x_m = genon_mum[,x]
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
  }))
  
  ## Determine the segregation of the groups
  group <- list()
  group$BI <- which(config == 1)
  group$PI <- which(config %in% c(2,3))
  group$MI <- which(config %in% c(4,5))

  obj <- list(genon=genon, depth_Ref=depth_Ref,depth_Alt=depth_Alt,chrom=chrom,pos=pos,indID=indID,SNP_Names=SNP_Names,config=config,
              group=group, nInd=nInd, nSnps=nSnps, filtering=filtering)
  
  ## Add the filtering at some stage
  obj$filtering <- NULL
    
  ## Compute the 2-point recombination fraction estimates
  obj <- rf_2pt(obj, nClust=nClust)
  
  ## return an object of the type we want
  return(obj)
}

dataG <- readRA(genofile,"Tassel",mum,dad)

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
  rf.PI <- foreach(snp1 = iter(1:nSnps_PI), .packages='GUSMap', .combine=comb) %do% {
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
  rf.MI <- foreach(snp1 = iter(1:nSnps_MI), .packages='GUSMap', .combine=comb) %dopar% {
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
  rf.BI <- foreach(snp1 = iter(1:nSnps_BI), .packages='GUSMap', .combine=comb) %dopar% {
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
      rf.ind <- switch(which.min(c(rf.est1[[2]], rf.est2[[2]], rf.est4[[2]]) ),
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
  rf.PI.BI <- foreach(snp.ps = iter(1:nSnps_PI), .packages='GUSMap', .combine=comb) %dopar% {
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
  rf.MI.BI <- foreach(snp.ms = iter(1:nSnps_MI), .packages='GUSMap', .combine=comb) %dopar% {
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



createGroups <- function(obj, LOD=4){
  
  ## Initialize the LGs
  LGs <- list()
  
  ## Start with the maternal side
  unmapped <- sort(c(obj$group$BI,obj$group$MI))

  ## Run the process of computing the linkage groups
  finish = FALSE
  while(!finish){
    newLG <- unmapped[sort(which(obj$LOD[unmapped,unmapped]==min(obj$LOD[unmapped,unmapped],na.rm=T),arr.ind=T)[1,])]
    unmapped <- unmapped[-which(unmapped%in%newLG)]
    compLOD <- obj$LOD[newLG,unmapped]
    minLOD <- min(compLOD)
    while(minLOD < LODthres ){
      newSNP <- matrix(which(compLOD == minLOD,arr.ind=T),ncol=2)[1,]
      if(newSNP[1]<length(newLG))
        newLG <- c(newLG[1:newSNP[1]],unmapped[newSNP[2]],newLG[(newSNP[1]+1):length(newLG)])
      else
        newLG <- c(newLG,unmapped[newSNP[2]])
      unmapped <- unmapped[-newSNP[2]]
      compLOD <- obj$LOD[newLG,unmapped]
      minLOD <- min(compLOD)
    }
    finish <- !any(obj$LOD[unmapped,unmapped] < LODthres)
    LGs <- c(LGs,list(newLG))
  }
  
  obj$LGs <- LGs
  return(obj)
}


mergeGroups <- function(obj, LOD=3){
  
}




