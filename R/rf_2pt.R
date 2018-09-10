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
#' FS method: Compute 2-point recombination fraction estimates and LOD scores
#' 
#' Method for estimating 2-point recombination fraction and associated LOD scores.
#' 
#' In this function, the recombination fraction for all pairs of SNPs is computed. Estimation 
#' is performed by optimimizing the hidden Markov model (HMM) likelihood given by 
#' \insertCite{bilton2018genetics1;textual}{GUSMap} for two markers, where the parental phase 
#' that is used is the one that maximizes the likelihood, given the segregation type that has been 
#' inferred. The \code{err} specifies whether sequencing errors should
#' also be estimated when computing the 2-point recombination fraction estimates. LOD scores associated 
#' with each pair of SNPs are also computed.
#' 
#' Parallelization is used to speed up the estimation process using the \code{\link[foreach]{foreach}}
#' function. The argument \code{nClust} specifies the number of cores to use in the parallelization. 
#' Note: It is important that the number of cores you specify is not more than what is available on 
#' your comupter (otherwise bad things can happen!).
#' 
#' Currently, only 2-point recombination fractions for a single full-sib family can be computed. However,
#' there are plans to extend this function to multiple full-sib families in the future.
#'  
#' Note: This function can take a while, espically when there are a large number of SNPs. 
#' 
#' @usage
#' FSobj$rf_2pt(nClust = 2, err = FALSE)
#' 
#' @param nClust An integer value for the number of cores to use in the parallelization of 
#' computing the 2-point recombination fraction estimates. 
#' @param err Locial value. If \code{TRUE}, a sequencing error parameter is also estimated as 
#' part of the 2-point recombination fraction estimates. If \code{FALSE}, then the sequencing 
#' error parameter is fixed at zero. 
#' 
#' @name $rf_2pt
#' @seealso \code{\link{FS}}
#' @author Timothy P. Bilton
#' @references 
#' \insertRef{bilton2018genetics1}{GUSMap}
#' @examples
#' ## simulate some sequencing data
#' set.seed(6745)
#' config <- list(list(sample(c(1,2,4), size=30, replace=T)))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50)
#' 
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt(nClust=1)
#' @aliases $rf_2pt


## function needed for foreach loop
comb <- function(...){
  mapply('rbind',...,SIMPLIFY=FALSE)
}

### Function for computing the pairwise RF in a single full-sib family.
rf_2pt_single <- function(ref, alt, config, config_infer, group, group_infer, nClust, nInd, err=FALSE){
  
  nSnps = as.integer(2)
  noFam = as.integer(1)
  init_ep = ifelse(err, 0.001, 0)
  ## if estimating sequencing error. Make sure that the initial value is nonzero
  if(err & init_ep < 0.0001)
    init_ep <- 0.0001
  
  if(length(c(unlist(group), unlist(group_infer))) == 0)
    stop("There are no SNPs in the data set.")
  
  ## Calcuate the number of SNPs
  indx_BI <- c(group$BI, group_infer$BI)
  nSnps_BI <- length(indx_BI)
  indx_MI <- c(group$MI, group_infer$SI)
  nSnps_MI <- length(indx_MI)
  indx_PI <- c(group$PI)
  nSnps_PI = length(indx_PI)
  
  ## set up the config vector
  if(!is.null(config_infer))
    config[which(!is.na(config_infer))] <- config_infer[which(!is.na(config_infer))]
  ## Check that there are no missing configs in the data
  if(any(is.na(config[c(indx_MI , indx_BI , indx_PI)])))
    stop("There are some missing segregation types in the data.")
  
  ## Set up the Clusters
  cl <- makeCluster(nClust)
  registerDoSNOW(cl)

  cat("\nComputing 2-point recombination fraction estimates ...\n")
  cat("Paternal informative SNPs\n")
  ## Paternal informative SNPs
  rf.PI <- foreach(snp1 = iter(seq(length.out=nSnps_PI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_PI),simplify=F)
    for(snp2 in seq_len(snp1-1)){
      ind = indx_PI[c(snp1,snp2)]
      ref2 <- list(ref[,ind])
      alt2 <- list(alt[,ind])
      #bcoef_mat <- list(choose(ref2[[1]]+alt2[[1]],ref2[[1]]))
      #Kab       <- list(bcoef_mat[[1]]*(1/2)^(ref2[[1]]+alt2[[1]]))
      ## compute the rf's and LOD
      rf.est1 <- rf_est_FS_2pt(0.2, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(5,5) + 2*(config[ind]==3))), seqErr=err)
      rf.est2 <- rf_est_FS_2pt(0.2, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(5,6) + 2*(config[ind]==3))), seqErr=err)
      rf.ind1 <- switch(which.max(c(rf.est1$loglik,rf.est2$loglik)), TRUE, FALSE)
      if(rf.ind1){
        rf[[1]][snp2] <- rf.est1$rf
        rf[[2]][snp2] <- rf.est1$LOD
      }
      else{
        rf[[1]][snp2] <- rf.est2$rf
        rf[[2]][snp2] <- rf.est2$LOD
      }
      # rf.est1 <- optim(GUSbase::logit2(0.2), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(5,5) + 2*(config[ind]==3))), seqErr=err)
      # rf.est2 <- optim(GUSbase::logit2(0.2), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(5,6) + 2*(config[ind]==3))), seqErr=err)
      # rf.ind1 <- switch(which.max(c(rf.est1$loglik,rf.est2$loglik)), TRUE, FALSE)
      # if(rf.ind1){
      #   rf[[1]][snp2] <- rf.est1$rf
      #   if(rf.est1$rf > 0.499)
      #     rf[[2]][snp2] <- 0
      #   else
      #     rf[[2]][snp2] <- rf.est1$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                           nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                           OPGP=list(as.integer(c(5,5) + 2*(config[ind]==3))),
      #                                           extra=rf.est2$ep)
      # } else{
      #   rf[[1]][snp2] <- rf.est2$rf
      #   if(rf.est2$rf > 0.499)
      #     rf[[2]][snp2] <- 0
      #   else
      #     rf[[2]][snp2] <- rf.est2$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(5,6) + 2*(config[ind]==3))),
      #                                             extra=rf.est2$ep)
      # }
    }
    return(rf)
  }
  for(i in 1:2){
    rf.PI[[i]][upper.tri(rf.PI[[i]])] <- t(rf.PI[[i]])[upper.tri(rf.PI[[i]])]
  }
  
  cat("Maternal informative SNPs\n")
  ## Maternal informative SNPs
  rf.MI <- foreach(snp1 = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_MI),simplify=F)
    for(snp2 in seq_len(snp1-1)){
      ind = indx_MI[c(snp1,snp2)]
      ref2 <- list(ref[,ind])
      alt2 <- list(alt[,ind])
      #bcoef_mat <- list(choose(ref2[[1]]+alt2[[1]],ref2[[1]]))
      #Kab       <- list(bcoef_mat[[1]]*(1/2)^(ref2[[1]]+alt2[[1]]))
      ## compute the rf's and LOD
      rf.est1 <- rf_est_FS_2pt(0.2, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(9,9) + 2*(config[ind]==5))), seqErr=err)
      rf.est2 <- rf_est_FS_2pt(0.2, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(9,10) + 2*(config[ind]==5))), seqErr=err)
      rf.ind1 <- switch(which.max(c(rf.est1$loglik,rf.est2$loglik)), TRUE, FALSE)
      if(rf.ind1){
        rf[[1]][snp2] <- rf.est1$rf
        rf[[2]][snp2] <- rf.est1$LOD
      }
      else{
        rf[[1]][snp2] <- rf.est2$rf
        rf[[2]][snp2] <- rf.est2$LOD
      }
      #rf.est1 <- optim(GUSbase::logit2(0.2), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(9,9) + 2*(config[ind]==5))), seqErr=F)
      # rf.est2 <- optim(GUSbase::logit2(0.2), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(9,10) + 2*(config[ind]==5))), seqErr=F)
      # rf.ind1 <- switch(which.max(c(rf.est1$loglik,rf.est2$loglik)), TRUE, FALSE)
      # if(rf.ind1){
      #   rf[[1]][snp2] <- rf.est1$rf
      #   if(rf.est1$rf > 0.499)
      #     rf[[2]][snp2] <- 0
      #   else
      #     rf[[2]][snp2] <- rf.est1$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(9,9) + 2*(config[ind]==5))),
      #                                             extra=rf.est1$ep)
      # } else{
      #   rf[[1]][snp2] <- rf.est2$rf
      #   if(rf.est2$rf > 0.499)
      #     rf[[2]][snp2] <- 0
      #   else
      #     rf[[2]][snp2] <- rf.est2$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(9,10) + 2*(config[ind]==5))),
      #                                             extra=rf.est2$ep)
      # }
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
      ref2 <- list(ref[,ind])
      alt2 <- list(alt[,ind])
      #bcoef_mat <- list(choose(ref2[[1]]+alt2[[1]],ref2[[1]]))
      #Kab       <- list(bcoef_mat[[1]]*(1/2)^(ref2[[1]]+alt2[[1]]))
      ## compute the rf's and LOD
      # phase 1
      temp1 <- rf_est_FS_2pt(0.1, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(1,1))), seqErr=err)
      temp2 <- rf_est_FS_2pt(0.4, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(1,1))), seqErr=err)
      # temp1 <- optim(GUSbase::logit2(0.1), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(1,1))), seqErr=F)
      # temp2 <- optim(GUSbase::logit2(0.4), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                OPGP=list(as.integer(c(1,1))), seqErr=F)
      rf.est1 <- switch(which.max(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
      # phase 2
      temp1 <- rf_est_FS_2pt(0.1, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(1,2))), seqErr=err)
      temp2 <- rf_est_FS_2pt(0.4, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(1,2))), seqErr=err)
      # temp1 <- optim(GUSbase::logit2(0.1), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                OPGP=list(as.integer(c(1,2))), seqErr=F)
      # temp2 <- optim(GUSbase::logit2(0.4), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                OPGP=list(as.integer(c(1,2))), seqErr=F)
      rf.est2 <- switch(which.max(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
      # phase 3
      temp1 <- rf_est_FS_2pt(0.1, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(1,4))), seqErr=err)
      temp2 <- rf_est_FS_2pt(0.4, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(1,4))), seqErr=err)
      # temp1 <- optim(GUSbase::logit2(0.1), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                OPGP=list(as.integer(c(1,4))), seqErr=F)
      # temp2 <- optim(GUSbase::logit2(0.4), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                OPGP=list(as.integer(c(1,4))), seqErr=F)
      rf.est4 <- switch(which.max(c(temp1$loglik,temp2$loglik)),temp1,temp2) 
      ## work out which is best
      rf.ind <- switch(which.max(c(rf.est1$loglik,rf.est2$loglik,rf.est4$loglik)), 1, 2, 3)
      if(rf.ind == 1){
        rf[[1]][snp2] <- rf.est1$rf
        rf[[2]][snp2] <- rf.est1$LOD
      }
      else if(rf.ind == 2){
        rf[[1]][snp2] <- rf.est2$rf
        rf[[2]][snp2] <- rf.est2$LOD
      }
      else {
        rf[[1]][snp2] <- rf.est4$rf
        rf[[2]][snp2] <- rf.est4$LOD
      }
      #   if(rf.est1$rf > 0.499)
      #     rf[[2]][snp2] <- 0
      #   else
      #     rf[[2]][snp2] <- rf.est1$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(1,1))),
      #                                             extra=rf.est1$ep)
      # } else if(rf.ind == 2){
      #   rf[[1]][snp2] <- rf.est2$rf
      #   if(rf.est2$rf > 0.499)
      #     rf[[2]][snp2] <- 0
      #   else
      #     rf[[2]][snp2] <- rf.est2$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(1,2))),
      #                                             extra=rf.est2$ep)
      # } else{
      #   rf[[1]][snp2] <- rf.est4$rf
      #   if(rf.est4$rf > 0.499)
      #     rf[[2]][snp2] <- 0
      #   else
      #     rf[[2]][snp2] <- rf.est4$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(1,4))),
      #                                             extra=rf.est4$ep)
      # }
    }
    return(rf)
  }
  for(i in 1:2){
    rf.BI[[i]][upper.tri(rf.BI[[i]])] <- t(rf.BI[[i]])[upper.tri(rf.BI[[i]])]
  }
  
  cat("Paternal informative vs Both informative\n")
  ## Paternal and Informative SNPs
  rf.PI.BI <- foreach(snp.ps = iter(seq(length.out=nSnps_PI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_BI),simplify=F)
    for(snp.bi in 1:nSnps_BI){
      ind <- c(indx_PI[snp.ps],indx_BI[snp.bi])
      ref2 <- list(ref[,ind])
      alt2 <- list(alt[,ind])
      #bcoef_mat <- list(choose(ref2[[1]]+alt2[[1]],ref2[[1]]))
      #Kab       <- list(bcoef_mat[[1]]*(1/2)^(ref2[[1]]+alt2[[1]]))
      ## compute the rf's and LOD
      rf.est1 <- rf_est_FS_2pt(0.2, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(5,1) + 2*c(config[ind[[1]]]==3,0))), seqErr=err)
      rf.est2 <- rf_est_FS_2pt(0.2, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(5,2) + 2*c(config[ind[[1]]]==3,0))), seqErr=err)
      # rf.est1 <- optim(GUSbase::logit2(0.2), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(5,1) + 2*c(config[ind[[1]]]==3,0))), seqErr=F)
      # rf.est2 <- optim(GUSbase::logit2(0.2), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(5,2) + 2*c(config[ind[[1]]]==3,0))), seqErr=F)
      rf.ind1 <- switch(which.max(c(rf.est1$loglik,rf.est2$loglik)), TRUE, FALSE)
      if(rf.ind1){
        rf[[1]][snp.bi] <- rf.est1$rf
        rf[[2]][snp.bi] <- rf.est1$LOD
      }
      else{
        rf[[1]][snp.bi] <- rf.est2$rf
        rf[[2]][snp.bi] <- rf.est2$LOD
      }
      #   if(rf.est1$rf > 0.499)
      #     rf[[2]][snp.bi] <- 0
      #   else
      #     rf[[2]][snp.bi] <- rf.est1$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(5,1) + 2*c(config[ind[[1]]]==3,0))),
      #                                             extra=rf.est1$ep)
      # } else{
      #   rf[[1]][snp.bi] <- rf.est2$rf
      #   if(rf.est2$rf > 0.499)
      #     rf[[2]][snp.bi] <- 0
      #   else
      #     rf[[2]][snp.bi] <-  rf.est2$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(5,2) + 2*c(config[ind[[1]]]==3,0))),
      #                                             extra=rf.est2$ep)
      # }
    }
    return(rf)
  }
  
  cat("Maternal informative vs Both informative\n")
  ## Maternal and Informative SNPs
  rf.MI.BI <- foreach(snp.mi = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_BI),simplify=F)
    for(snp.bi in 1:nSnps_BI){
      ind <- c(indx_MI[snp.mi],indx_BI[snp.bi])
      ref2 <- list(ref[,ind])
      alt2 <- list(alt[,ind])
      #bcoef_mat <- list(choose(ref2[[1]]+alt2[[1]],ref2[[1]]))
      #Kab       <- list(bcoef_mat[[1]]*(1/2)^(ref2[[1]]+alt2[[1]]))
      ## compute the rf's and LOD
      rf.est1 <- rf_est_FS_2pt(0.2, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(9,1) + 2*c(config[ind[[1]]]==5,0))), seqErr=err)
      rf.est2 <- rf_est_FS_2pt(0.2, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(9,3) + 2*c(config[ind[[1]]]==5,0))), seqErr=err)
      # rf.est1 <- optim(GUSbase::logit2(0.2), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(9,1) + 2*c(config[ind[[1]]]==5,0))), seqErr=F)
      # rf.est2 <- optim(GUSbase::logit2(0.2), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(9,3) + 2*c(config[ind[[1]]]==5,0))), seqErr=F)
      rf.ind1 <- switch(which.max(c(rf.est1$loglik,rf.est2$loglik)), TRUE, FALSE)
      if(rf.ind1){
        rf[[1]][snp.bi] <- rf.est1$rf
        rf[[2]][snp.bi] <- rf.est1$LOD
      } else{
        rf[[1]][snp.bi] <- rf.est2$rf
        rf[[2]][snp.bi] <- rf.est2$LOD
      }
      #   if(rf.est1$rf > 0.499)
      #     rf[[2]][snp.bi] <- 0
      #   else
      #     rf[[2]][snp.bi] <- rf.est1$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(9,1) + 2*c(config[ind[[1]]]==5,0))),
      #                                             extra=rf.est1$ep)
      # } else{
      #   rf[[1]][snp.bi] <- rf.est2$rf
      #   if(rf.est2$rf > 0.499)
      #     rf[[2]][snp.bi] <- 0
      #   else
      #     rf[[2]][snp.bi] <- rf.est2$loglik +
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(9,3) + 2*c(config[ind[[1]]]==5,0))),
      #                                             extra=rf.est2$ep)
      # }
    }
    return(rf)
  }
  
  ## For the non-informative computations
  ## Really done so that we can check that there is no miss identification of the group
  cat("Maternal informative vs Paternal informative\n")
  rf.MI.PI <- foreach(snp.mi = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_PI),simplify=F)
    for(snp.pi in 1:nSnps_PI){
      ind <- c(indx_MI[snp.mi],indx_PI[snp.pi])
      ref2 <- list(ref[,ind])
      alt2 <- list(alt[,ind])
      #bcoef_mat <- list(choose(ref2[[1]]+alt2[[1]],ref2[[1]]))
      #Kab       <- list(bcoef_mat[[1]]*(1/2)^(ref2[[1]]+alt2[[1]]))
      ## compute the rf's and LOD
      rf.est1 <- rf_est_FS_2pt(0.2, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(9,9) + 2*c(config[ind] %in% c(3,5)))), seqErr=err)
      rf.est2 <- rf_est_FS_2pt(0.2, ep=init_ep, ref=ref2,alt=alt2, OPGP=list(as.integer(c(9,10) + 2*c(config[ind] %in% c(3,5)))), seqErr=err)
      # rf.est1 <- optim(GUSbase::logit2(0.2), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(9,9) + 2*c(config[ind] %in% c(3,5)))), seqErr=F)
      # rf.est2 <- optim(GUSbase::logit2(0.2), fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err, method="BFGS",
      #                  ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                  nInd=nInd,nSnps=nSnps,noFam=noFam,
      #                  OPGP=list(as.integer(c(9,10) + 2*(config[ind] %in% c(3,5)))), seqErr=F)
      rf.ind1 <- switch(which.max(c(rf.est1$loglik,rf.est2$loglik)), TRUE, FALSE)
      if(rf.ind1){
        rf[[1]][snp.pi] <- rf.est1$rf
        rf[[2]][snp.pi] <- rf.est1$LOD
      } else{
        rf[[1]][snp.pi] <- rf.est2$rf
        rf[[2]][snp.pi] <- rf.est2$LOD
      }
      #   if(rf.est1$rf > 0.499)
      #     rf[[2]][snp.pi] <- 0
      #   else
      #     rf[[2]][snp.pi] <- rf.est1$value + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(9,9) + 2*c(config[ind] %in% c(3,5)))),
      #                                             extra=rf.est$ep)
      # } else{
      #   rf[[1]][snp.pi] <- rf.est2$rf
      #   if(rf.est2$rf > 0.499)
      #     rf[[2]][snp.pi] <- 0
      #   else
      #     rf[[2]][snp.pi] <- rf.est2$loglik + 
      #                         ll_fs_mp_scaled_err(1000,ref=ref2,alt=alt2,bcoef_mat=bcoef_mat,Kab=Kab,
      #                                             nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
      #                                             OPGP=list(as.integer(c(5,2) + 2*c(config[ind] %in% c(3,5))))))
      # }
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
  LOD.mat[which(LOD.mat < 0)] <- 0  # some negative LOD scores from rounding
  return(list(rf = rf.mat,LOD = LOD.mat))
}



### Function for computing the pairwise RF in a multiple-sib family.
rf_2pt_multi <- function(ref, alt, config, group, nClust, noFam, init_r = 0.25){
  
  if(length(unlist(group)) == 0)
    stop("There are no SNPs in the data set.")
  
  ## Calcuate the number of SNPs
  indx_BI <- c(group$BI)
  nSnps_BI <- length(indx_BI)
  indx_MI <- c(group$MI)
  nSnps_MI <- length(indx_MI)
  indx_PI <- c(group$PI)
  nSnps_PI = length(indx_PI)
  
  ## Set up the Clusters
  cl <- makeCluster(nClust)
  registerDoSNOW(cl)
  
  cat("\nComputing 2-point recombination fraction estimates ...\n")
  cat("Paternal informative SNPs\n")
  ## Paternal informative SNPs
  rf.PI <- foreach(snp1 = iter(seq(length.out=nSnps_PI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_PI),simplify=F)
    for(snp2 in seq_len(snp1-1)){
      ind = indx_PI[c(snp1,snp2)]
      configFam <- matrix(unlist(lapply(config, function(z) z[ind])),ncol=2, nrow=noFam,byrow=T)
      wFam <- which(apply(configFam,1, function(z) all(z %in% c(2,3))))
      OPGP <- vector(mode="list", length=noFam)
      ## Work out the phase
      for(fam in wFam){
        if(all(configFam[fam,] %in% c(2,3))){
          OPGP1 <- c(5,5) + 2*(configFam[fam,]==3)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(5,6) + 2*(configFam[fam,]==3)
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), OPGP1, OPGP2)
        }
        else if(all(configFam[fam,]==1)){
          OPGP1 <- c(1,1)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                      OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1,2)
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                      OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP3 <- c(1,4)
          rf.est3 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP3), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik,rf.est3$loglik)), OPGP1, OPGP2, OPGP3)
        }
        else if(any(configFam[fam,] == 1) & any(configFam[fam,] %in% c(2,3))){
          OPGP1 <- c(1,1) + 4*(configFam[fam,]==2) + 6*(configFam[fam,]==3)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1,2) + 4*(configFam[fam,]==2) + 6*(configFam[fam,]==3)
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), OPGP1, OPGP2)
        }
      }
      rf.est <- GUSMap:::rf_est_FS(init_r=init_r,ref=lapply(ref[wFam], function(x) x[,ind]),
                                   alt=lapply(alt[wFam], function(x) x[,ind]),
                                   OPGP=OPGP[wFam], epsilon=NULL, noFam=length(wFam), nThreads=1)
      rf[[1]][snp2] <- rf.est$rf
      rf[[2]][snp2] <- rf.est$LOD
    }
    return(rf)
  }
  for(i in 1:2){
    rf.PI[[i]][upper.tri(rf.PI[[i]])] <- t(rf.PI[[i]])[upper.tri(rf.PI[[i]])]
  }
  
  cat("Maternal informative SNPs\n")
  ## Maternal informative SNPs
  rf.MI <- foreach(snp1 = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_MI),simplify=F)
    for(snp2 in seq_len(snp1-1)){
      ind = indx_MI[c(snp1,snp2)]
      configFam <- matrix(unlist(lapply(config, function(z) z[ind])),ncol=2, nrow=noFam,byrow=T)
      wFam <- which(apply(configFam,1, function(z) all(z %in% c(4,5))))
      OPGP <- vector(mode="list", length=noFam)
      ## Work out the phase
      for(fam in wFam){
        if(all(configFam[fam,] %in% c(4,5))){
          OPGP1 <- c(9,9) + 2*(configFam[fam,]==5)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(9,10) + 2*(configFam[fam,]==5)
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), OPGP1, OPGP2)
        }
        else if(all(configFam[fam,]==1)){
          OPGP1 <- c(1,1)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1,2)
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP3 <- c(1,4)
          rf.est3 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP3), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik,rf.est3$loglik)), OPGP1, OPGP2, OPGP3)
        }
        else if(any(configFam[fam,] == 1) & any(configFam[fam,] %in% c(4,5))){
          OPGP1 <- c(1,1) + 8*(configFam[fam,]==4) + 10*(configFam[fam,]==5)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1 + 8*(configFam[fam,1]==4) + 10*(configFam[fam,1]==5),3 + 7*(configFam[fam,2]==4) + 9*(configFam[fam,2]==5))
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), OPGP1, OPGP2)
        }
      }
      rf.est <- GUSMap:::rf_est_FS(init_r=init_r,ref=lapply(ref[wFam], function(x) x[,ind]),
                                   alt=lapply(alt[wFam], function(x) x[,ind]),
                                   OPGP=OPGP[wFam], epsilon=NULL, noFam=length(wFam), nThreads=1)
      rf[[1]][snp2] <- rf.est$rf
      rf[[2]][snp2] <- rf.est$LOD
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
      configFam <- matrix(unlist(lapply(config, function(z) z[ind])),ncol=2, nrow=noFam,byrow=T)
      wFam <- which(apply(configFam,1, function(z) (all(z %in% 1)) | (z[1] %in% c(2:5) &  z[2] == 1) | (z[2] %in% c(2:5) &  z[1] == 1)))
      OPGP <- vector(mode="list", length=noFam)
      ## Work out the phase
      for(fam in wFam){
        if(all(configFam[fam,]==1)){
          OPGP1 <- c(1,1)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1,2)
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP3 <- c(1,4)
          rf.est3 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP3), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik,rf.est3$loglik)), OPGP1, OPGP2, OPGP3)
        }
        else if( any(configFam[fam,]==1) & any(configFam[fam,]%in%c(2,3)) ){
          OPGP1 <- c(1,1) + 4*(configFam[fam,]==2) + 6*(configFam[fam,]==3)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1,2) + 4*(configFam[fam,]==2) + 6*(configFam[fam,]==3)
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), OPGP1, OPGP2)
        }
        else if( any(configFam[fam,]==1) & any(configFam[fam,]%in%c(4,5)) ){
          OPGP1 <- c(1,1) + 8*(configFam[fam,]==4) + 10*(configFam[fam,]==5)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1 + 8*(configFam[fam,1]==4) + 10*(configFam[fam,1]==5),3 + 7*(configFam[fam,2]==4) + 9*(configFam[fam,2]==5))
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), OPGP1, OPGP2)
        }
      }
      rf.est <- GUSMap:::rf_est_FS(init_r=init_r,ref=lapply(ref[wFam], function(x) x[,ind]),
                                   alt=lapply(alt[wFam], function(x) x[,ind]),
                                   OPGP=OPGP[wFam], epsilon=NULL, noFam=length(wFam), nThreads=1)
      rf[[1]][snp2] <- rf.est$rf
      rf[[2]][snp2] <- rf.est$LOD
    }
    return(rf)
  }
  for(i in 1:2){
    rf.BI[[i]][upper.tri(rf.BI[[i]])] <- t(rf.BI[[i]])[upper.tri(rf.BI[[i]])]
  }
  
  cat("Paternal informative vs Both informative\n")
  ## Paternal and Informative SNPs
  rf.PI.BI <- foreach(snp.pi = iter(seq(length.out=nSnps_PI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_BI),simplify=F)
    for(snp.bi in 1:nSnps_BI){
      ind <- c(indx_PI[snp.pi],indx_BI[snp.bi])
      configFam <- matrix(unlist(lapply(config, function(z) z[ind])),ncol=2, nrow=noFam,byrow=T)
      wFam <- which(apply(configFam,1, function(z) z[1] %in% c(2,3) & z[2] == 1))
      OPGP <- vector(mode="list", length=noFam)
      ## Work out the phase
      for(fam in wFam){
        if(all(configFam[fam,]==1)){
          OPGP1 <- c(1,1)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1,2)
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP3 <- c(1,4)
          rf.est3 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP3), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik,rf.est3$loglik)), OPGP1, OPGP2, OPGP3)
        }
        else if(any(configFam[fam,] == 1) & any(configFam[fam,] %in% c(2,3))){
          OPGP1 <- c(1,1) + 4*(configFam[fam,]==2) + 6*(configFam[fam,]==3)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1,2) + 4*(configFam[fam,]==2) + 6*(configFam[fam,]==3)
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), OPGP1, OPGP2)
        }
      }
      rf.est <- GUSMap:::rf_est_FS(init_r=init_r,ref=lapply(ref[wFam], function(x) x[,ind]),
                                   alt=lapply(alt[wFam], function(x) x[,ind]),
                                   OPGP=OPGP[wFam], epsilon=NULL, noFam=length(wFam), nThreads=1)
      rf[[1]][snp.bi] <- rf.est$rf
      rf[[2]][snp.bi] <- rf.est$LOD
    }
    return(rf)
  }
  
  cat("Maternal informative vs Both informative\n")
  ## Maternal and Informative SNPs
  rf.MI.BI <- foreach(snp.mi = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_BI),simplify=F)
    for(snp.bi in 1:nSnps_BI){
      ind <- c(indx_MI[snp.mi],indx_BI[snp.bi])
      configFam <- matrix(unlist(lapply(config, function(z) z[ind])),ncol=2, nrow=noFam,byrow=T)
      wFam <- which(apply(configFam,1, function(z) z[1] %in% c(4,5) & z[2] == 1))
      OPGP <- vector(mode="list", length=noFam)
      ## Work out the phase
      for(fam in wFam){
        if(all(configFam[fam,]==1)){
          OPGP1 <- c(1,1)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1,2)
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP3 <- c(1,4)
          rf.est3 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP3), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik,rf.est3$loglik)), OPGP1, OPGP2, OPGP3)
        }
        else if(any(configFam[fam,] == 1) & any(configFam[fam,] %in% c(4,5))){
          OPGP1 <- c(1,1) + 8*(configFam[fam,]==4) + 10*(configFam[fam,]==5)
          rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
          OPGP2 <- c(1 + 8*(configFam[fam,1]==4) + 10*(configFam[fam,1]==5),3 + 7*(configFam[fam,2]==4) + 9*(configFam[fam,2]==5))
          rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                        OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
          OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), OPGP1, OPGP2)
        }
      }
      rf.est <- GUSMap:::rf_est_FS(init_r=init_r,ref=lapply(ref[wFam], function(x) x[,ind]),
                                    alt=lapply(alt[wFam], function(x) x[,ind]),
                                    OPGP=OPGP[wFam], epsilon=NULL, noFam=length(wFam), nThreads=1)
      rf[[1]][snp.bi] <- rf.est$rf
      rf[[2]][snp.bi] <- rf.est$LOD
    }
    return(rf)
  }
  
  ## For the non-informative computations
  ## Really done so that we can check that there is no miss identification of the group
  cat("Maternal informative vs Paternal informative\n")
  rf.MI.PI <- foreach(snp.mi = iter(seq(length.out=nSnps_MI)), .combine=comb) %dopar% {
    rf <- replicate(2,numeric(nSnps_PI),simplify=F)
    for(snp.pi in 1:nSnps_PI){
      ind = c(indx_MI[snp.mi],indx_PI[snp.pi])
      configFam <- matrix(unlist(lapply(config, function(z) z[ind])),ncol=2, nrow=noFam,byrow=T)
      wFam <- which(apply(configFam,1, function(z) z[1] %in% c(4,5) & z[2] %in% c(2,3)))
      OPGP <- vector(mode="list", length=noFam)
      ## Work out the phase
      for(fam in wFam){
        OPGP1 <- c(5,5) + 2*(configFam[fam,] %in% c(3,5))
        rf.est1 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                      OPGP=list(OPGP1), epsilon=NULL, nThreads=1)
        OPGP2 <- c(5,6) + 2*(configFam[fam,] %in% c(3,5))
        rf.est2 <- GUSMap:::rf_est_FS(init_r=init_r,ref=list(ref[[fam]][,ind]),alt=list(alt[[fam]][,ind]),
                                      OPGP=list(OPGP2), epsilon=NULL, nThreads=1)
        OPGP[[fam]] <- switch(which.min(c(rf.est1$loglik,rf.est2$loglik)), OPGP1, OPGP2)
      }
      rf.est <- GUSMap:::rf_est_FS(init_r=init_r,ref=lapply(ref[wFam], function(x) x[,ind]),
                                   alt=lapply(alt[wFam], function(x) x[,ind]),
                                   OPGP=OPGP[wFam], epsilon=NULL, noFam=length(wFam), nThreads=1)
      rf[[1]][snp.pi] <- rf.est$rf
      rf[[2]][snp.pi] <- rf.est$LOD
    }
    return(rf)
  }
  
  #rf.MI.PI <- replicate(2, matrix(NA, nrow=nSnps_MI, ncol=nSnps_PI),simplify=FALSE)
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




rf_est_FS_2pt <- function(init_r=0.01, ep=0.001, ref, alt, OPGP,
                          seqErr=T, trace=F, noFam=as.integer(1), method = "optim", nThreads = 1){
  
  nInd <- lapply(ref,nrow)  # number of individuals
  nSnps <- as.integer(2)   # number of SNPs
  
  if(method=="optim"){
    
    ## Compute the K matrix for heterozygous genotypes
    bcoef_mat <- Kab <- vector(mode="list", length=noFam)
    for(fam in 1:noFam){
      bcoef_mat[[fam]] <- choose(ref[[fam]]+alt[[fam]],ref[[fam]])
      Kab[[fam]] <- bcoef_mat[[fam]]*(1/2)^(ref[[fam]]+alt[[fam]])
    }
    
    if(seqErr)
      para <- c(GUSbase::logit2(init_r),GUSbase::logit(ep))
    else
      para <- GUSbase::logit2(init_r)

    ## Find MLE
    optim.MLE <- optim(para, fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err,
                       method="BFGS", 
                       ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                       nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam,
                       seqErr=seqErr,extra=ep,nThreads=nThreads)
    ep_hat <- ifelse(seqErr,GUSbase::inv.logit(optim.MLE$par[2]),ep)
    ## find the LOD score
    LOD <- -(optim.MLE$value - ll_fs_mp_scaled_err(1000,ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                                                nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
                                                OPGP=OPGP, extra=ep_hat, nThreads=nThreads))
    # return the results
    return(list(rf=GUSbase::inv.logit2(optim.MLE$par[1]), 
                ep=ep_hat, loglik=-optim.MLE$value, LOD=LOD))
  } else{ # EM algorithm approach
    ## Set up the parameter values
    EM.arg = c(500,1e-15)
    
    ## convert the data into the right format:
    OPGPmat = do.call(what = "rbind",OPGP)
    ref_mat = do.call(what = "rbind",ref)
    alt_mat = do.call(what = "rbind",alt)
    
    EMout <- .Call("EM_HMM", rep(init_r,2), ep, ref_mat, alt_mat, OPGPmat,
                   noFam, unlist(nInd), nSnps, FALSE, seqErr, EM.arg, as.integer(0), nThreads=nThreads)
    
    EMout[[3]] = EMout[[3]] + sum(log(choose(ref_mat+alt_mat,ref_mat)))
    ## find the LOD score
    LOD <- EMout[[3]] + ll_fs_mp_scaled_err(1000,ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                                                    nInd=nInd,nSnps=nSnps,noFam=noFam,seqErr=F,
                                                    OPGP=OPGP, extra=EMout[[2]], nThreads=nThreads)
    return(list(rf=EMout[[1]][1], ep=EMout[[2]],
                loglik=EMout[[3]], LOD=LOD))
  }
}

