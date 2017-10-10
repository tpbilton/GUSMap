#### Functions for converting phase information to OPGPs
#### Author: Timothy P. Bilton
#### Date: 03/04/17
#### Edited: 19/06/17


## Function for deteriming the parental phase of a full-sib family
infer_OPGP_FS <- function(depth_Ref, depth_Alt, config, epsilon=NULL, delta=NULL, ...){
  
  if(!is.matrix(depth_Ref) || !is.matrix(depth_Alt))
    stop("The read counts inputs are not matrix objects")
  if( (!is.null(epsilon) & !is.numeric(epsilon)) || (!is.null(delta) & !is.numeric(delta)) )
    stop("Starting values for the error parameters needs to be a single numeric value in the interval (0,1) or a NULL object")
  nSnps <- ncol(depth_Ref); nInd <- nrow(depth_Ref)
  
  ## solve the likelihood for when phase is not known
  MLEs <- rf_est_FS_UP(depth_Ref, depth_Alt, config, epsilon=epsilon, delta=delta, ...)
  
  parHap <- matrix("A",nrow=4,ncol=nSnps)
  ## Determine whether the phase between adjacent paternal (or maternal) informative SNPs is in coupling
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

## Function for converting parental haplotypes to OPGPs
parHapToOPGP <- function(parHap, major="A", minor="B"){
  return (apply(parHap,2,function(x){
    if(sum(x==major)==2){
      return((x[1]==minor)+(x[3]==minor)*2+1)
    } else if (x[1]==x[2] & x[3]==x[4]){
      return(9)
    } else if (all(x[3:4]==major)){
      return((x[1]==minor)+5)
    } else if (all(x[1:2]==major)){
      return((x[3]==minor)+7)
    } 
  }) 
)}



