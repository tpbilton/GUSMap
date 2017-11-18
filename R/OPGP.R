#### Functions for converting phase information to OPGPs
#### Author: Timothy P. Bilton
#### Date: 03/04/17
#### Edited: 19/06/17


## Function for deteriming the parental phase of a full-sib family
infer_OPGP_FS <- function(depth_Ref, depth_Alt, config, epsilon=0.001, ...){
  
  if(!is.matrix(depth_Ref) || !is.matrix(depth_Alt))
    stop("The read counts inputs are not matrix objects")
  if( (!is.null(epsilon) & !is.numeric(epsilon)) )
    stop("Starting values for the error parameter needs to be a single numeric value in the interval (0,1) or a NULL object")
  if( !is.numeric(config) || !is.vector(config) || any(!(config == round(config))) || any(config < 1) || any(config > 9) )
    stop("Segregation information needs to be an integer vector equal to the number of SNPs with entires from 1 to 9")
  
  nSnps <- ncol(depth_Ref); nInd <- nrow(depth_Ref)
  
  ## Index the informative loci
  Isnps <- which(!(config %in% 6:9))

  ## solve the likelihood for when phase is not known
  MLEs <- rf_est_FS_UP(depth_Ref[,Isnps], depth_Alt[,Isnps], config[Isnps], epsilon=epsilon, ...)
  
  parHap <- matrix("A",nrow=4,ncol=nSnps)
  parHap[cbind(c(rep(3,sum(config %in% c(3,7,9))),rep(4,sum(config %in% c(3,7,9)))),which(config %in% c(3,7,9)))] <- "B"
  parHap[cbind(c(rep(1,sum(config %in% c(5,8,9))),rep(2,sum(config %in% c(5,8,9)))),which(config %in% c(5,8,9)))] <- "B"
  ## Determine whether the phase between adjacent paternal (or maternal) informative SNPs is in coupling
  coupling_pat <- MLEs[[1]] < 0.5
  coupling_mat <- MLEs[[2]] < 0.5
  
  # Work out the indices of the r.f. parameters of each sex
  ps_snps <- which(config %in% c(1,2,3))
  ms_snps <- which(config %in% c(1,4,5))
  npar <- c(ifelse(length(ps_snps)==0,0,length(ps_snps)-1),
            ifelse(length(ms_snps)==0,0,length(ms_snps)-1))

  if(any(npar==0)){ ## essentially a backcross approach
    ## If BC on maternal side
    if(npar[1] == 0){
      s_mat = numeric(npar[2] + 1)
      s_mat[1] <- 2
      for(i in 1:npar[2])
        s_mat[i+1] <- ifelse(coupling_mat[i],s_mat[i],s_mat[i]%%2+1)
      parHap[cbind(s_mat+2,ms_snps)] <- "B"
    }
    ## If BC on paternal side
    if(npar[2] == 0){
      s_pat = numeric(npar[1] + 1)
      s_pat[1] <- 2
      for(i in 1:npar[1])
        s_pat[i+1] <- ifelse(coupling_pat[i],s_pat[i],s_pat[i]%%2+1)
      parHap[cbind(s_pat,ps_snps)] <- "B"
    }
  }
  else{ # Some intercross SNPs are present
    s_pat <- numeric(npar[1] + 1)
    s_pat[1] <- 2
    
    s_mat <- numeric(npar[2] + 1)
    s_mat[1] <- 2
    
    for(i in 1:npar[1]){
      s_pat[i+1] <- ifelse(coupling_pat[i],s_pat[i],s_pat[i]%%2+1)
    }
    for(i in 1:npar[2]){
      s_mat[i+1] <- ifelse(coupling_mat[i],s_mat[i],s_mat[i]%%2+1)
    }
    
    parHap[cbind(s_pat,ps_snps)] <- "B"
    parHap[cbind(s_mat+2,ms_snps)] <- "B"
    
  }
  
  return(parHapToOPGP(parHap))
}

## Function for converting parental haplotypes to OPGPs
parHapToOPGP <- function(parHap, major="A", minor="B"){
  return (apply(parHap,2,function(x){
    if (all(x==major)){
      return(13) 
    } else if (x[1]==x[2] & x[3]==x[4] & x[1]==major & x[3]==minor){
      return(14)
    } else if (x[1]==x[2] & x[3]==x[4] & x[3]==major & x[1]==minor){
      return(15)
    } else if (all(x==minor)){
      return(16)
    } else if(sum(x==major)==2){
      return((x[1]==minor)+(x[3]==minor)*2+1)
    } else if (all(x[3:4]==major)){
      return((x[1]==minor)+5)
    } else if (all(x[1:2]==major)){
      return((x[3]==minor)+9)
    } else if (all(x[3:4]==minor)){
      return((x[1]==minor)+7)
    } else if (all(x[1:2]==minor)){
      return((x[3]==minor)+11)
    } 
  }) 
)}



