### Wrapper functions required for calling C code.
### Author: Timothy Bilton
### Date: 18/01/17
### Edited: 19/06/17

#### Wrapper function for likelihoods of the GBS HMM when
#### Likelihood function is written in C in the file 'likelihoods.c'

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

## r.f.'s are sex-specific and constrained to the range [0,1] (for unphased data)
ll_fs_up_ss_scaled <- function(logit_r,genon,depth,config,nInd,nSnps,ps,ms,npar){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit(logit_r[1:npar[1]])
  r[ms,2] <- inv.logit(logit_r[npar[1]+1:npar[2]])
  .Call("ll_fs_up_ss_scaled_c",r,genon,depth,config,nInd,nSnps)
}
