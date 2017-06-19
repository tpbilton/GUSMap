#### Functions for converting phase information to OPGPs
#### Author: Timothy P. Bilton
#### Date: 03/04/17


## Function for converting parental haplotypes to OPGPs
parHapToOPGP <- function(parHap){
  return (apply(parHap,2,function(x){
    if(sum(x=='A')==2){
      return((x[1]=='B')+(x[3]=='B')*2+1)
    } else if (sum(x=='A')==4){
      return(0)
    } else if (all(x[3:4]=='A')){
      return((x[1]=='B')+5)
    } else if (all(x[1:2]=='A')){
      return((x[3]=='B')+7)
    }
  }) 
)}

## Function for converting phase information in the form of inheritance vectors
## into to OPGPs
phaseToOPGP <- function(phase, nSnps, swap=FALSE){
  parHap <- matrix("A",nrow=4,ncol=nSnps)
  if(swap){
    mat <- substr(phase,1,1)
    pat <- substr(phase,2,2)
  }
  else{
    pat <- substr(phase,1,1)
    mat <- substr(phase,2,2)
  }
  
  ## Sort out any duplicates
  if(any(pat=='d')|any(mat=='u')){
    pat[which(pat == 'd')] <- pat[which(pat == 'd')-1]
    mat[which(mat == 'd')] <- pat[which(mat == 'd')-1]
  }
  
  parHap[cbind(suppressWarnings(as.numeric(pat))+1,1:nSnps)] <- "B"
  parHap[cbind(suppressWarnings(as.numeric(mat))+3,1:nSnps)] <- "B"
  
  ## make sure the baseline phase is correct
  if(parHap[1,which(apply(parHap[1:2,],2,function(x) !(all(x=='A'))))[1]] == 'B')
    parHap[1:2,] <- parHap[2:1,]
  if(parHap[3,which(apply(parHap[3:4,],2,function(x) !(all(x=='A'))))[1]] == 'B')
    parHap[3:4,] <- parHap[4:3,]
  
  ## Convert the parental haplotypes into phase numbers
  OPGP <- parHapToOPGP(parHap)
}

## Function for deteriming the parental phase of a full-sib family
phaseDet_FS <- function(genon,depth,config, scaled=F, ...){
  nSnps <- ncol(genon); nInd <- nrow(genon)
  parHap <- matrix("A",nrow=4,ncol=nSnps)
  
  ## solve the likelihood for when phase is not known
  MLEs <- rf_est_FS_UP(genon, depth, config, scaled=scaled, ...)
  
  ## Determine whether the phase between adjacent paternal (or maternal) segregating SNPs is in coupling
  coupling_pat <- MLEs[[1]] < 0.5
  coupling_mat <- MLEs[[2]] < 0.5
  
  # Work out the indices of the r.f. parameters of each sex
  ps_snps <- which(config %in% c(1,2))
  ms_snps <- which(config %in% c(1,3))
  npar <- c(length(ps_snps)-1,length(ms_snps)-1)
  
  s_pat <- numeric(npar[[1]]+1)
  s_pat[1] <- 2
  
  s_mat <- numeric(npar[[2]]+1)
  s_mat[1] <- 2
  
  for(i in 1:npar[[1]]){
    s_pat[i+1] <- ifelse(coupling_pat[i],s_pat[i],s_pat[i]%%2+1)
  }
  for(i in 1:npar[[2]]){
    s_mat[i+1] <- ifelse(coupling_mat[i],s_mat[i],s_mat[i]%%2+1)
  }
  
  parHap[cbind(s_pat,ps_snps)] <- "B"
  parHap[cbind(s_mat+2,ms_snps)] <- "B"
  
  return(parHapToOPGP(parHap))
}

