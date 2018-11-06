##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping
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
#' FS method: Compute ordered linkage maps
#' 
#' Method for inferring parental phase (e.g., ordered parental genotype pair (OPGP)) and
#' estimating recombination fractions in full-sib families.
#' 
#' This function infers the parental phase (or ordered parental genotype pair (OPGP)) and
#' estimates adjacent recombination fractions using the hidden Markov model (HMM) approach as
#' described in \insertCite{bilton2018genetics1;textual}{GUSMap}. 
#' 
#' The optimization of the likelihood for the HMM is performed using either the Expectation-Maximumization (EM) algorithm
#' (\code{method="EM"}) or using direct numeric optimization via the \code{\link{optim}} function (\code{method="optim"}).
#' The likelihood computations (and computation of derivatives if required) are scaled using 
#' forward and backward recursion to avoid overflow issues and are performed in C. These computations 
#' are also parallelization via the OpenMP package, where the argument \code{nThreads} specifies
#' how many threads to use. Be careful not to set \code{nThreads} to more than the number of threads available
#' on your computer (or bad things will happen). In addition, if the package is complied without OpenMP, then this 
#' parallelization has no effect and the likelihood is computed in serial.
#' 
#' If \code{mapped = TRUE}, then combined linkage groups must have been formed from the \code{\link{$addBIsnps}} function
#' first (and preferably ordered from the \code{\link{$orderLG}} function).
#' 
#' @usage 
#' FSobj$computeMap(chrom=NULL, init_r=0.01, ep=0.001, method="optim", sexSpec=FALSE, 
#'                     err=TRUE, mapped=TRUE, nThreads=1)
#' 
#' @param chrom A integer vector giving the indices of the chromosomes (or linkage groups) to be computed.
#' @param init_r A numeric value giving the initial values for the recombination fractions. Each 
#' recombination fraction parameter is set to the same initial value.
#' @param ep A numeric value giving the initial value for the sequencing error parameter.
#' @param method A character string specifying whether optimization should be performed using
#' direct maximization (\code{optim}) or via the Expectation-Maximum (EM) algorithm (\code{EM}).
#' @param sexSpec Logical value. If \code{TRUE}, sex-specific recombination fractions are
#' are estimated.
#' @param err Locical value. If \code{TRUE}, the sequencing error parameter is estimated. Otherwise, the
#' sequenicng error parameter is fixed to the value of the \code{ep} argument.
#' @param mapped Locial value. If \code{TRUE}, the maps are computed using the marker order given 
#' in the combined linkage groups. Otherwise, the maps are computed using the original marker order  
#' given by the genomic assembly.
#' @param nThreads An integer value giving the number of threads to use in computing the likelihood in parallel.
#' @name $computeMap
#' @author Timothy P. Bilton and Chris Scott
#' @seealso \code{\link{FS}}
#' @references
#' \insertRef{bilton2018genetics1}{GUSMap}
#' @examples
#' #### Case 1: Compute linkage map from linkage groups
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(list(sample(c(1,2,4), size=30, replace=TRUE)))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50, epsilon=0.005)
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt(nClust=1)
#' ## create and order linkage groups
#' F1data$createLG()
#' F1data$addBIsnps()
#' F1data$orderLG(ndim=5)
#' 
#' ## Compute the linkage map
#' F1data$computeMap()
#' 
#' #### Case 2: Compute map using original assembly order
#' F1data$computeMap(mapped = FALSE)
#' @aliases $computeMap

rf_est_FS <- function(init_r=0.01, ep=0.001, ref, alt, OPGP, noFam=as.integer(1),
                      sexSpec=FALSE, seqErr=TRUE, multiErr=FALSE, trace=FALSE,
                      method = "optim", nThreads=1, ...){
  
  ## Do some checks
  nInd <- lapply(ref,nrow)  # number of individuals
  nSnps <- ncol(ref[[1]])   # number of SNPs

  if(method=="optim"){
    # Arguments for the optim function
    optim.arg <- list(...)
    if(length(optim.arg) == 0)
      optim.arg <- list(maxit = 5000, reltol=1e-15)
    
    ## Compute the K matrix for heterozygous genotypes
    bcoef_mat <- Kab <- vector(mode="list", length=noFam)
    for(fam in 1:noFam){
      bcoef_mat[[fam]] <- choose(ref[[fam]]+alt[[fam]],ref[[fam]])
      Kab[[fam]] <- bcoef_mat[[fam]]*(1/2)^(ref[[fam]]+alt[[fam]])
    }
    
    ## If we want to estimate sex-specific r.f.'s
    if(sexSpec){
      
      stop("To be implemented")
      
      # Work out the indices of the r.f. parameters of each sex
      ps <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% 1:8)))))[-1] - 1
      ms <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% c(1:4,9:12))))))[-1] - 1
      npar <- c(length(ps),length(ms))
      
      # Determine the initial values
      if(length(init_r)==1) 
        para <- GUSbase::logit2(rep(init_r,sum(npar)))
      else if(length(init_r) != sum(npar)) 
        para <- GUSbase::logit2(rep(0.1,sum(npar)))
      else
        para <- init_r
      # sequencing error
      if(length(ep) != 1 & !is.null(ep))
        para <- c(para,GUSbase::logit(0.001))
      else if(!is.null(ep))
        para <- c(para,GUSbase::logit(ep))
      
      ## Are we estimating the error parameters?
      seqErr=!is.null(ep)
      
      ## Find MLE
      optim.MLE <- stats::optim(para,ll_fs_ss_mp_scaled_err,method="BFGS",control=optim.arg,
                         ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                         nInd=nInd,nSnps=nSnps,OPGP=OPGP,ps=ps,ms=ms,npar=npar,noFam=noFam,
                         seqErr=!is.null(ep))
    }
    else{
      # Determine the initial values
      if(length(init_r)==1) 
        para <- GUSbase::logit2(rep(init_r,nSnps-1))
      else if(length(init_r) != nSnps-1) 
        para <- GUSbase::logit2(rep(0.1,nSnps-1))
      else
        para <- init_r
      # sequencing error
      if(seqErr){
        if(multiErr) para <- c(para,GUSbase::logit(rep(ep, nSnps)))
        else para <- c(para,GUSbase::logit(ep))
      }
      
      ## Find MLE
      optim.MLE <- stats::optim(para, fn=ll_fs_mp_scaled_err, gr=score_fs_mp_scaled_err,
                         method="BFGS", control=optim.arg,
                         ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                         nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam,
                         seqErr=seqErr,extra=ep,nThreads=nThreads,multiErr=multiErr)
    }
    # Print out the output from the optim procedure (if specified)
    if(trace){
      print(optim.MLE)
    }
    # Work out what to return for the recombination fractions
    if(sexSpec)
      rfReturn <- list(rf_p=GUSbase::inv.logit2(optim.MLE$par[1:npar[1]]),
                       rf_m=GUSbase::inv.logit2(optim.MLE$par[npar[1]+1:npar[2]]))
    else
      rfReturn <- list(rf=GUSbase::inv.logit2(optim.MLE$par[1:(nSnps-1)]))
    # work out what to return for the sequencing errors
    if(seqErr){
      if(sexSpec){
        if(multiErr) epReturn <- list(GUSbase::inv.logit(optim.MLE$par[sum(npar) + 1:(2*nSnps)]))
        else         epReturn <- list(GUSbase::inv.logit(optim.MLE$par[sum(npar) + 1]))
      } else{
        if(multiErr) epReturn <- list(GUSbase::inv.logit(optim.MLE$par[nSnps:(2*nSnps-1)]))
        else         epReturn <- list(GUSbase::inv.logit(optim.MLE$par[nSnps]))
      }
    }
    else epReturn <- ep
    # Check for convergence
    if(trace & optim.MLE$convergence != 0)
      warning(paste0('Optimization failed to converge properly with error ',optim.MLE$convergence,'\n smallest MLE estimate is: ', round(min(optim.MLE$par),6)))
    # Return the results
    return(c(rfReturn, ep=epReturn, loglik=-optim.MLE$value))
  } else{ # EM algorithm approach
    if(multiErr)
      stop("Yet to be implemented")
    else{
      ## Set up the parameter values
      temp.arg <- list(...)
      if(!is.null(temp.arg$maxit) && is.numeric(temp.arg$maxit) && length(temp.arg$maxit) == 1) 
        EM.arg = c(temp.arg$maxit)
      else
        EM.arg = c(1000)
      if(!is.null(temp.arg$reltol) && is.numeric(temp.arg$reltol) && length(temp.arg$reltol) == 1){
        EM.arg = c(EM.arg,temp.arg$reltol)
      }
      else
        EM.arg = c(EM.arg,1e-20)
      
      # Determine the initial values
      if(length(init_r)==1)
        init_r <- rep(init_r,2*(nSnps-1))
      else if(length(init_r) != nSnps-1)
        init_r <- rep(0.1,2*(nSnps-1))
      
      if(sexSpec){
        ps <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% 1:8)))))[-1] - 1
        ms <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% c(1:4,9:12))))))[-1] - 1
        npar <- c(length(ps),length(ms))
        ss_rf <- logical(2*(nSnps-1))
        ss_rf[ps] <- TRUE
        ss_rf[ms + nSnps-1] <- TRUE
      }
      else ss_rf = 0;
      ## convert the data into the right format:
      OPGPmat = do.call(what = "rbind",OPGP)
      ref_mat = do.call(what = "rbind",ref)
      alt_mat = do.call(what = "rbind",alt)
      
      EMout <- .Call("EM_HMM", init_r, ep, ref_mat, alt_mat, OPGPmat,
                     noFam, unlist(nInd), nSnps, sexSpec, seqErr, EM.arg, as.integer(ss_rf), nThreads=nThreads)
      
      EMout[[3]] = EMout[[3]] + sum(log(choose(ref_mat+alt_mat,ref_mat)))
      
      if(sexSpec){
        return(list(rf_p=EMout[[1]][ps],rf_m=EMout[[1]][nSnps-1+ms],
                    ep=EMout[[2]],
                    loglik=EMout[[3]]))
      }
      else
        return(list(rf=EMout[[1]][1:(nSnps-1)], 
                    ep=EMout[[2]],
                    loglik=EMout[[3]]))
      
    }
  }
}


## recombination estimates for case where the phase is unkonwn.
## The r.f.'s are sex-specific and constrained to the range [0,1]
rf_est_FS_UP <- function(ref, alt, config, ep, method="optim", trace=F, nThreads=0, ...){
  
  ## Check imputs
  if( any( ref<0 | !is.finite(ref)) || any(!(ref == round(ref))))
    stop("The read count matrix for the reference allele is invalid")
  if( any( alt<0 | !is.finite(alt)) || any(!(alt == round(alt))))
    stop("The read count matrix for the alternate allele is invalid")
  if( !is.vector(config) || !is.numeric(config) || any(!(config %in% 1:9)) )
    stop("Invalid config vector. It must be a numeric vector with entries between 1 and 9.")
  if( !is.logical(trace) || is.na(trace) )
    trace = FALSE
  if(!(method %in% c("EM","optim")))
    stop("Specified optimization method is unknown. Please select one of 'EM' or 'optim'")
  
  nInd <- nrow(ref)  # number of individuals
  nSnps <- ncol(ref)  # number of SNPs
  if(!is.integer(ref))
    ref <- matrix(as.integer(ref), nrow=nInd, ncol=nSnps)
  if(!is.integer(alt))
    alt <- matrix(as.integer(alt), nrow=nInd, ncol=nSnps)
  if(!is.integer(config))
    config <- as.integer(config)
  
  if(method == "optim"){
    # Arguments for the optim function
    optim.arg <- list(...)
    if(length(optim.arg) == 0)
      optim.arg <- list(maxit = 1000, reltol=1e-10)
    
    # Work out the indices of the r.f. parameters of each sex
    ps <- which(config %in% c(1,2,3))[-1] - 1
    ms <- which(config %in% c(1,4,5))[-1] - 1
    npar <- c(length(ps),length(ms))
    
    ## Compute the K matrix for heterozygous genotypes
    bcoef_mat <- choose(ref+alt,ref)
    Kab <- bcoef_mat*(1/2)^(ref+alt)
    
    ## Are we estimating the error parameters?
    seqErr=!is.null(ep)
    
    para <- GUSbase::logit(rep(0.5,sum(npar)))
    # sequencing error
    if(length(ep) != 1 & !is.null(ep))
      para <- c(para,GUSbase::logit(0.01))
    else if(!is.null(ep))
      para <- c(para,GUSbase::logit(ep))
    
    if(nSnps > 2){
      ## Find MLE
      optim.MLE <- stats::optim(para,ll_fs_up_ss_scaled_err,method="BFGS",control=optim.arg,
                         ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                         nInd=nInd,nSnps=nSnps,config=config,ps=ps,ms=ms,npar=npar,
                         seqErr=seqErr,nThreads=nThreads)
      # Print out the output from the optim procedure (if specified)
      if(trace){
        print(optim.MLE)
      }
      
      # Check for convergence
      if(optim.MLE$convergence != 0)
        warning(paste0('Optimization failed to converge properly with error ',optim.MLE$convergence,'\n smallest MLE estimate is: ', round(min(optim.MLE$par),6)))
      
      # Return the MLEs
      return(list(rf_p=GUSbase::inv.logit(optim.MLE$par[1:npar[1]]),rf_m=GUSbase::inv.logit(optim.MLE$par[npar[1]+1:npar[2]]),
                  ep=GUSbase::inv.logit(optim.MLE$par[sum(npar)+1])))
    } 
    else if(nSnps == 2){
      ## If both SNPs are informative, need to use the Nelder-Mead to distinguish between the two sexes.
      if(all(config == 1)){
        optim.MLE <- stats::optim(para,ll_fs_up_ss_scaled_err,method="Nelder-Mead",control=optim.arg,
                           ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                           nInd=nInd,nSnps=nSnps,config=config,ps=ps,ms=ms,npar=npar,
                           seqErr=seqErr)
      }
      ## Otherwise, proceed as normal
      else{
        optim.MLE <- stats::optim(para,ll_fs_up_ss_scaled_err,method="BFGS",control=optim.arg,
                           ref=ref,alt=alt,bcoef_mat=bcoef_mat,Kab=Kab,
                           nInd=nInd,nSnps=nSnps,config=config,ps=ps,ms=ms,npar=npar,
                           seqErr=seqErr)
      }
      # Print out the output from the optim procedure (if specified)
      if(trace){
        print(optim.MLE)
      }
      
      # Check for convergence
      if(trace & optim.MLE$convergence != 0)
        warning(paste0('Optimization failed to converge properly with error ',optim.MLE$convergence,'\n smallest MLE estimate is: ', round(min(optim.MLE$par),6)))
      
      # Return the MLEs
      return(list(rf_p=inv.logit(optim.MLE$par[npar[1]]),rf_m=inv.logit(optim.MLE$par[npar[1]+npar[2]]),
                  ep=inv.logit(optim.MLE$par[sum(npar)+1])))
    }
  }
  else{ ## EM approach
    ## Set up the parameter values
    temp.arg <- list(...)
    if(!is.null(temp.arg$maxit) && is.numeric(temp.arg$maxit) && length(temp.arg$maxit) == 1) 
      EM.arg = c(temp.arg$maxit)
    else
      EM.arg = c(5000)
    if(!is.null(temp.arg$reltol) && is.numeric(temp.arg$reltol) && length(temp.arg$reltol) == 1){
      EM.arg = c(EM.arg,temp.arg$reltol)
    }
    else
      EM.arg = c(EM.arg,1e-5)
    
    ## work out which rf can be estimated
    ps <- which(config %in% c(1,2,3))[-1] - 1
    ms <- which(config %in% c(1,4,5))[-1] - 1 
    npar <- c(length(ps),length(ms))
    ss_rf <- logical(2*(nSnps-1))
    ss_rf[ps] <- TRUE
    ss_rf[ms + nSnps-1] <- TRUE
    
    ## Are we estimating the error parameters?
    seqErr=!is.null(ep)
    
    EMout <- .Call("EM_HMM_UP", rep(0.5,(nSnps-1)*2), ep, ref, alt, config,
                   as.integer(1), nInd, nSnps, seqErr, EM.arg, as.integer(ss_rf), nThreads=nThreads)
    return(list(rf_p=EMout[[1]][ps],rf_m=EMout[[1]][nSnps-1+ms],
                ep=EMout[[2]],
                loglik=EMout[[3]]))
  }
}



