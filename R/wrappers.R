##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2018 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################
### Wrapper functions required for calling C code.
### Author: Timothy Bilton
#' @useDynLib GUSMap

#### Likelihood function is written in C in the file 'likelihoods.c'
## r.f.'s are equal
ll_fs_mp_scaled_err <- function(para,ref,alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam,seqErr,extra=0,nThreads=1,multiErr=FALSE){
  ## untransform the parameters
  r <- GUSbase::inv.logit2(para[1:(nSnps-1)])
  if(seqErr){
    if(multiErr)  ep = GUSbase::inv.logit(para[nSnps:(2*nSnps-1)])
    else          ep = rep(GUSbase::inv.logit(para[nSnps]), length.out = nSnps)
  } else ep = rep(extra, length.out = nSnps)
  ## define likelihood
  llval = 0
  for(fam in 1:noFam)
    llval = llval + .Call("ll_fs_scaled_err_c",r,ep,ref[[fam]], alt[[fam]], bcoef_mat[[fam]], Kab[[fam]],
                          OPGP[[fam]],nInd[[fam]],nSnps,nThreads)
  return(llval)
}

## r.f.'s are sex-specific
ll_fs_ss_mp_scaled_err <- function(para,ref,alt,bcoef_mat,Kab,OPGP,nInd,nSnps,ps,ms,npar,noFam,seqErr){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- GUSbase::inv.logit2(para[1:npar[1]])
  r[ms,2] <- GUSbase::inv.logit2(para[npar[1]+1:npar[2]])
  if(seqErr)
    ep = GUSbase::inv.logit(para[sum(npar)+1])
  else
    ep = 0
  ## define likelihood
  llval = 0
  # define the density values for the emission probs
  Kaa <- Kbb <- vector(mode = "list", length=noFam)
  for(fam in 1:noFam){
    Kaa[[fam]] <- bcoef_mat[[fam]]*(1-ep)^ref[[fam]]*ep^alt[[fam]]
    Kbb[[fam]] <- bcoef_mat[[fam]]*(1-ep)^alt[[fam]]*ep^ref[[fam]]
  }
  for(fam in 1:noFam)
    llval = llval + .Call("ll_fs_ss_scaled_err_c",r,Kaa[[fam]],Kab[[fam]],Kbb[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
  return(llval)
}

## r.f.'s are sex-specific and constrained to the range [0,1] (for unphased data)
ll_fs_up_ss_scaled_err <- function(para,ref,alt,bcoef_mat,Kab,config,nInd,nSnps,ps,ms,npar,seqErr,nThreads=1){
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  r[ps,1] <- GUSbase::inv.logit(para[1:npar[1]])
  r[ms,2] <- GUSbase::inv.logit(para[npar[1]+1:npar[2]])
  if(seqErr)
    ep = GUSbase::inv.logit(para[sum(npar)+1])
  else
    ep = 0
  ## define likelihood
  llval = 0
  .Call("ll_fs_up_ss_scaled_err_c",r,bcoef_mat,ep,ref,alt,Kab,config,nInd,nSnps,nThreads)
}

#### Score functions written in C in the file 'score.c'
## r.f.'s are equal
score_fs_mp_scaled_err <- function(para,ref,alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam,seqErr=TRUE,multiErr=FALSE,extra=0,nThreads=1){
  ## untransform the parameters
  r <- GUSbase::inv.logit2(para[1:(nSnps-1)])
  ## Compute score vector
  llval=0
  if(seqErr){
    # Sequencing error parameter for each SNP
    if(multiErr){
      score = numeric(2*nSnps-1)
      ep <- GUSbase::inv.logit(para[nSnps-1 + 1:nSnps])
      for(fam in 1:noFam){
        temp = .Call("score_fs_scaled_multi_err_c",r,ep,ref[[fam]],alt[[fam]],
                     bcoef_mat[[fam]],Kab[[fam]],OPGP[[fam]],nInd[[fam]],nSnps,nThreads)
        llval = llval + temp[[1]]
        score = score + temp[[2]]
      }
    }
    # A single sequencing error parameter for all the SNPs
    else{
      score = numeric(nSnps)
      ep <- GUSbase::inv.logit(para[nSnps])
      for(fam in 1:noFam){
        temp = .Call("score_fs_scaled_err_c",r,ep,ref[[fam]],alt[[fam]],
                     bcoef_mat[[fam]],Kab[[fam]],OPGP[[fam]],nInd[[fam]],nSnps,nThreads)
        llval = llval + temp[[1]]
        score = score + temp[[2]]
      }
    }
  } else{  # No sequencing error parameter
    score = numeric(nSnps - 1)
    ep = extra
    for(fam in 1:noFam){
      temp = .Call("score_fs_scaled_c",r,ep,ref[[fam]], alt[[fam]], bcoef_mat[[fam]],
                   Kab[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
      llval = llval + temp[[1]]
      score = score + temp[[2]]
    }
  }
  assign(".score",-score, envir = parent.frame(3))
  return(llval)
}

score_fs_mp_ss_scaled_err <- function(para,ref,alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam,ps,ms,npar,seqErr=TRUE,multiErr=FALSE,extra=0,nThreads=1){
  ## untransform the parameters
  r <- matrix(0,ncol=2,nrow=nSnps-1)
  if(npar[1] > 0) r[ps,1] <- GUSbase::inv.logit2(para[1:npar[1]])
  if(npar[2] > 0) r[ms,2] <- GUSbase::inv.logit2(para[npar[1]+1:npar[2]])
  if(seqErr){
    llval = 0
    if(!multiErr){
      score = numeric(2*nSnps-1)
      ep <- GUSbase::inv.logit(para[length(para)])
      for(fam in 1:noFam){
        psf <- which(OPGP[[fam]] %in% 1:8)[-1] - 1
        msf <- which(OPGP[[fam]] %in% c(1:4,9:12))[-1] - 1
        ss_rf <- logical(2*(nSnps-1))
        ss_rf[psf] <- TRUE
        ss_rf[msf + nSnps-1] <- TRUE
        temp <- .Call("score_fs_ss_scaled_err_c",r,ep,ref[[fam]],alt[[fam]],
                      bcoef_mat[[fam]],Kab[[fam]],OPGP[[fam]],nInd[[fam]],nSnps,nThreads,as.integer(ss_rf))
        score = score + temp[[2]]
        llval = llval + temp[[1]]
      }
      ss_rf <- logical(2*nSnps-1)
      ss_rf[ps] <- TRUE
      ss_rf[ms + nSnps-1] <- TRUE
      ss_rf[2*nSnps-1] <- TRUE
      assign(".score",-score[which(ss_rf)], envir = parent.frame(3))
    } else{
      score = numeric(3*nSnps-2)
      ep <- GUSbase::inv.logit(para[(sum(npar)+1):length(para)])
      for(fam in 1:noFam){
        psf <- which(OPGP[[fam]] %in% 1:8)[-1] - 1
        msf <- which(OPGP[[fam]] %in% c(1:4,9:12))[-1] - 1
        ss_rf <- logical(2*(nSnps-1))
        ss_rf[psf] <- TRUE
        ss_rf[msf + nSnps-1] <- TRUE
        temp <- .Call("score_fs_ss_scaled_multi_err_c",r,ep,ref[[fam]],alt[[fam]],
                      bcoef_mat[[fam]],Kab[[fam]],OPGP[[fam]],nInd[[fam]],nSnps,nThreads,as.integer(ss_rf))
        score = score + temp[[2]]
        llval = llval + temp[[1]]
      }
      ss_rf <- logical(3*(nSnps-2))
      ss_rf[ps] <- TRUE
      ss_rf[ms + nSnps-1] <- TRUE
      ss_rf[(2*nSnps-1):(3*nSnps-2)] <- TRUE
      assign(".score",-score[which(ss_rf)], envir = parent.frame(3))
    }
  }
  return(llval)
}

## Function for extracting the score vector already computed
score_extract <- function(...){
  return(get(".score", envir = parent.frame(3)))
}
  

######### For computing the likelihood and derivatives at the same time:
# HMM_fs_mp_scaled_err <- function(para,ref,alt,bcoef_mat,Kab,OPGP,nInd,nSnps,noFam,seqErr){
#   ## untransform the parameters
#   r <- GUSbase:::inv.logit2(para[1:(nSnps-1)])
#   if(seqErr)
#     ep = GUSbase:::inv.logit(para[nSnps])
#   else
#     ep = 0
#   ## define likelihood
#   llval = 0
#   # define the density values for the emission probs
#   Kaa <- Kbb <- vector(mode = "list", length=noFam)
#   for(fam in 1:noFam){
#     Kaa[[fam]] <- bcoef_mat[[fam]]*(1-ep)^ref[[fam]]*ep^alt[[fam]]
#     Kbb[[fam]] <- bcoef_mat[[fam]]*(1-ep)^alt[[fam]]*ep^ref[[fam]]
#   }
#   for(fam in 1:noFam)
#     llval = llval + .Call("ll_fs_scaled_err_c",r,Kaa[[fam]],Kab[[fam]],Kbb[[fam]],OPGP[[fam]],nInd[[fam]],nSnps)
#   return(llval)
# }

