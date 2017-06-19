### Wrapper functions required for calling C code.
### Author: Timothy Bilton
### Date: 18/01/17
### Edited: 1/04/17

#### Wrapper function for likelihoods of the GBS HMM when
#### Likelihood function is written in C in the file 'GBS_LM_HMM_ll.c'
#  - logit2_r: Vector of transformed recombination fractions values 
#  - genon: The matrix of genotype calls. 0=BB, 1=AB, 2=AA, NA=missing
#  - depth: The matrix of depth calls. Positive for non-missing genotype calls.
#  - OPGP: Vector of OPGP's (that are assumed to be known or predetermined)
#  - nInd: Number of individuals
#  - nSnps: number of SNPs
#  - ps: Vector of indices for the paternal r.f. parameter to use in the transition matrix
#  - ms: Vector of indices for the maternal r.f. parameter to use in the transition matrix
#  - npar: Vector of number of r.f. parameters in each parent


## r.f.'s are equal
ll_fs <- function(logit2_r,genon,depth,OPGP,nInd,nSnps){
  r <- inv.logit2(logit2_r)
  .Call("ll_fs_c",r,genon,depth,OPGP,nInd,nSnps)
}

## r.f.'s are sex-specific
ll_fs_ss <- function(logit_r,genon,depth,OPGP,nInd,nSnps,ps,ms,npar){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit2(logit_r[1:npar[1]])
  r[ms,2] <- inv.logit2(logit_r[npar[1]+1:npar[2]])
  .Call("ll_fs_ss_c",r,genon,depth,OPGP,nInd,nSnps)
}

## r.f.'s are sex-specific and constrained to the range [0,1] (for unphased data)
ll_fs_up_ss <- function(logit_r,genon,depth,config,nInd,nSnps,ps,ms,npar){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit(logit_r[1:npar[1]])
  r[ms,2] <- inv.logit(logit_r[npar[1]+1:npar[2]])
  .Call("ll_fs_up_ss_c",r,genon,depth,config,nInd,nSnps)
}

#### For scaled likelihoods
## r.f.'s are equal
ll_fs_scaled <- function(logit2_r,genon,depth,OPGP,nInd,nSnps){
  r <- inv.logit2(logit2_r)
  .Call("ll_fs_scaled_c",r,genon,depth,OPGP,nInd,nSnps)
}

## r.f.'s are sex-specific
ll_fs_ss_scaled <- function(logit_r,genon,depth,OPGP,nInd,nSnps,ps,ms,npar){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit2(logit_r[1:npar[1]])
  r[ms,2] <- inv.logit2(logit_r[npar[1]+1:npar[2]])
  .Call("ll_fs_ss_scaled_c",r,genon,depth,OPGP,nInd,nSnps)
}

## r.f.'s are sex-specific and constrained to the range [0,1] (for unphased data)
ll_fs_up_ss_scaled <- function(logit_r,genon,depth,config,nInd,nSnps,ps,ms,npar){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit(logit_r[1:npar[1]])
  r[ms,2] <- inv.logit(logit_r[npar[1]+1:npar[2]])
  .Call("ll_fs_up_ss_scaled_c",r,genon,depth,config,nInd,nSnps)
}


##### likelihood functions for when there is more than 1 family

## r.f.'s are equal
ll_fs_mp <- function(logit2_r,genon,depth,OPGP,nInd,nSnps,noFam){
  r <- inv.logit2(logit2_r)
  llval = 0
  for(fam in 1:noFam)
    llval <- llVal + .Call("ll_fs_c",r,genon[[fam]],depth[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
  return(llval)
}

## r.f.'s are sex-specific
ll_fs_ss_mp <- function(logit_r,genon,depth,OPGP,nInd,nSnps,ps,ms,npar,noFam){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit2(logit_r[1:npar[1]])
  r[ms,2] <- inv.logit2(logit_r[npar[1]+1:npar[2]])
  llval = 0
  for(fam in 1:noFam)
    llval = llval + .Call("ll_fs_ss_c",r,genon[[fam]],depth[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
  return(llval)
}

## r.f.'s are equal
ll_fs_mp_scaled <- function(logit2_r,genon,depth,OPGP,nInd,nSnps,noFam){
  r <- inv.logit2(logit2_r)  
  llval = 0
  for(fam in 1:noFam)
    llval = llval + .Call("ll_fs_scaled_c",r,genon[[fam]],depth[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
  return(llval)
}

## r.f.'s are sex-specific
ll_fs_ss_mp_scaled <- function(logit_r,genon,depth,OPGP,nInd,nSnps,ps,ms,npar,noFam){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit2(logit_r[1:npar[1]])
  r[ms,2] <- inv.logit2(logit_r[npar[1]+1:npar[2]])
  llval = 0
  for(fam in 1:noFam)
    llval = llval + .Call("ll_fs_ss_scaled_c",r,genon[[fam]],depth[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
  return(llval)
}
