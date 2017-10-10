### Wrapper functions required for calling C code.
### Author: Timothy Bilton
### Date: 18/01/17
### Edited: 19/06/17

#### Wrapper function for likelihoods of the GBS HMM when
#### Likelihood function is written in C in the file 'likelihoods.c'

## r.f.'s are equal
ll_fs_mp_scaled_err <- function(para,depth_Ref,depth_Alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam,seqErr,allelicErr){
  ## untransform the parameters
  r <- inv.logit2(para[1:(nSnps-1)])
  if(seqErr)
    ep = inv.logit(para[nSnps])
  else
    ep = 0
  if(allelicErr)
    delta = as.numeric(inv.logit(para[length(para)]))
  else
    delta = 0
  ## define likelihood
  llval = 0
  # define the density values for the emission probs
  Kaa <- Kbb <- vector(mode = "list", length=noFam)
  for(fam in 1:noFam){
    Kaa[[fam]] <- bcoef_mat[[fam]]*(1-ep)^depth_Ref[[fam]]*ep^depth_Alt[[fam]]
    Kbb[[fam]] <- bcoef_mat[[fam]]*(1-ep)^depth_Alt[[fam]]*ep^depth_Ref[[fam]]
  }
  for(fam in 1:noFam)
    llval = llval + .Call("ll_fs_scaled_err_c",r,delta,Kaa[[fam]],Kab[[fam]],Kbb[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
  return(llval)
}

## r.f.'s are sex-specific
ll_fs_ss_mp_scaled_err <- function(para,depth_Ref,depth_Alt,bcoef_mat,Kab,OPGP,nInd,nSnps,ps,ms,npar,noFam,seqErr,allelicErr){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit2(para[1:npar[1]])
  r[ms,2] <- inv.logit2(para[npar[1]+1:npar[2]])
  if(seqErr)
    ep = inv.logit(para[nSnps])
  else
    ep = 0
  if(allelicErr)
    delta = as.numeric(inv.logit(para[length(para)]))
  else
    delta = 0
  ## define likelihood
  llval = 0
  # define the density values for the emission probs
  Kaa <- Kbb <- vector(mode = "list", length=noFam)
  for(fam in 1:noFam){
    Kaa[[fam]] <- bcoef_mat[[fam]]*(1-ep)^depth_Ref[[fam]]*ep^depth_Alt[[fam]]
    Kbb[[fam]] <- bcoef_mat[[fam]]*(1-ep)^depth_Alt[[fam]]*ep^depth_Ref[[fam]]
  }
  for(fam in 1:noFam)
    llval = llval + .Call("ll_fs_ss_scaled_err_c",r,delta,Kaa[[fam]],Kab[[fam]],Kbb[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
  return(llval)
}

## r.f.'s are sex-specific and constrained to the range [0,1] (for unphased data)
ll_fs_up_ss_scaled_err <- function(para,depth_Ref,depth_Alt,bcoef_mat,Kab,config,nInd,nSnps,ps,ms,npar,seqErr,allelicErr){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- inv.logit(para[1:npar[1]])
  r[ms,2] <- inv.logit(para[npar[1]+1:npar[2]])
  if(seqErr)
    ep = inv.logit(para[nSnps])
  else
    ep = 0
  if(allelicErr)
    delta = as.numeric(inv.logit(para[length(para)]))
  else
    delta = 0
  ## define likelihood
  llval = 0
  # define the density values for the emission probs
  Kaa <- bcoef_mat*(1-ep)^depth_Ref*ep^depth_Alt
  Kbb <- bcoef_mat*(1-ep)^depth_Alt*ep^depth_Ref
  .Call("ll_fs_up_ss_scaled_err_c",r,0,Kaa,Kab,Kbb,config,nInd,nSnps)
}
