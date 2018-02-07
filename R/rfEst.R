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
### R Script of functions for computing the recombination fraction estimate 
### for linkage analysis using sequencing data.
### Author: Timothy Bilton
### Date: 06/02/18

## Function for computing the recombination fraction when the parental phase is known
#' Estimation of adjacent recombination fractions in full-sib families.
#' 
#' Estimate the recombination fractions based on the hidden Markov model (HMM)
#' for low coverage sequencing data in full-sib families.
#' 
#' The \code{depth_Ref} and \code{depth_Alt} matrices must have the rows
#' representing the individuals and columns representing the SNPs. The entries
#' of these two matrices must be a finite non-negative integer number. The
#' number of families specified must match the number of genon (genotype calls)
#' and depth matrices given. For the vector of starting values, \code{init_r},
#' a single value can be given which sets all the starting values to be equal.
#' 
#' The likelihood calculations are scaled using forward recursion to avoid
#' overflow issues and so can be implemented on a large numbers of loci. If the
#' OPGP vector is unknown for the families, then it can be inferred using the function
#' \code{\link{infer_OPGP_FS}}.
#' 
#' Two different optimzation procedures are available, which are the EM algorithm and optim.
#' To control the parameters to these procedures, addition arguments can be passed to the function.
#' The arguments which have an effect are dependent on the optimization procedure.
#' \itemize{
#' \item EM: Only two arguments currently have an effect. 'reltol' specifies 
#' the maximum difference between the likelihood value of successive iterations
#' before the algorithm terminates. 'maxit' specifies the maximum number of iterations
#' used in the algorithm
#' \item optim: The extra arguments are passed directly to optim. Those see what 
#' arguments are valid, visit the help page fro optim using '?optim'.
#' }
#' 
#' @param init_r Vector of starting values for the recombination fractions
#' @param epsilon Numeric value of the starting value for the sequencing error
#' parameter. Default is NULL which means that the sequencing parameter is
#' fixed at zero.
#' @param depth_Ref List object with each element containing a matrix of allele
#' counts for the reference allele for each family.
#' @param depth_Alt List object with each element containing a matrix of allele
#' counts for the alternate allele for each family.
#' @param OPGP List object with each element containing a numeric vector of
#' ordered parental genotype pairs (OPGPs) for each family. See Bilton (2017)
#' for a classification of each OPGP.
#' @param sexSpec Logical value. If TRUE, sex specific recombination fractions
#' are estimated.
#' @param trace Logical value. If TRUE, output from \code{optim()} are printed.
#' @param noFam Numeric value. Specifies the number of full-sib families used
#' to estimate the recombination fractions.
#' @param method A character string specifying the optimzation procedure to be used.
#' @param \ldots Additional arguments passed to the optimizer procedure. See details for more information.
#' @return Function returns a list object. If non sex-specific recombination
#' fractions are specified, the list contains;
#' \itemize{
#' \item rf: Vector of recombination fraction estimates.
#' \item epsilon: Estimate of the sequencing error parameter.
#' \item loglik: The log-likelihood value at the maximum likelihood estimates.
#' }
#' else, if sex-specific recombination fractions are
#' desired, the list contains;
#' \itemize{
#' \item rf_p: Vector of the paternal recombination fraction estimates.
#' \item rf_m: Vector of the maternal recombination fraction estimates.
#' \item epsilon: Estimate of the sequencing error parameter.
#' \item loglik: The log-likelihood value at the maximum likelihood estimates.
#' }
#' @author Timothy P. Bilton
#' @seealso \code{\link{infer_OPGP_FS}}
#' @references Bilton, T.P., Schofield, M.R., Black, M.A., Chagné, D., Wilcox,
#' P.L., Dodds K.G. (2017). Accounting for errors in low coverage high-throughput
#' sequencing data when constructing genetic maps using biparental outcrossed
#' populations. Unpublished Manuscript.
#' @examples
#' 
#' ### Case 1: Single family
#' ## simulate full sib family
#' config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
#' F1data <- simFS(0.01, config=config, nInd=50, meanDepth=5)
#' 
#' ## Determine the parental phase
#' OPGP <- infer_OPGP_FS(F1data$depth_Ref, F1data$depth_Alt, config)
#' 
#' ## Estimate the recombination fractions
#' rf_est_FS(depth_Ref = list(F1data$depth_Ref), depth_Alt = list(F1data$depth_Alt), OPGP = list(OPGP), noFam = 1)
#' ## To change the optimzation parameters 
#' ## Max number of iterations for the EM algorithm set at 100
#' rf_est_FS(depth_Ref = list(F1data$depth_Ref), depth_Alt = list(F1data$depth_Alt), OPGP = list(OPGP),
#'   noFam = 1, maxit=100)
#' ## The algorithm will dtop when the difference between the likelihood at sucessive iterations is less
#' ## than 0.00001
#' rf_est_FS(depth_Ref = list(F1data$depth_Ref), depth_Alt = list(F1data$depth_Alt), OPGP = list(OPGP),
#'   noFam = 1, reltol=1e-5)
#' 
#' ########################
#' ### Case 2: Two families
#' config_1 <- c(1,1,3,1,4,4,1,6,2,4,6,1)
#' config_2 <- c(2,4,1,6,2,1,1,2,5,6,1,1)
#' 
#' Fam1 <- simFS(0.01, config=config_1, nInd=50, meanDepth=5)
#' Fam2 <- simFS(0.01, config=config_2, nInd=50, meanDepth=5)
#' 
#' ## Determine the parental phase in each family
#' OPGP_1 <- infer_OPGP_FS(Fam1$depth_Ref, Fam1$depth_Alt, config_1)
#' OPGP_2 <- infer_OPGP_FS(Fam2$depth_Ref, Fam2$depth_Alt, config_2)
#' 
#' ## Estimate the recombination fractions
#' rf_est_FS(depth_Ref = list(Fam1$depth_Ref,Fam2$depth_Ref),
#'           depth_Alt = list(Fam1$depth_Alt,Fam2$depth_Alt), OPGP = list(OPGP_1,OPGP_2), noFam = 2)
#'           
#' 
#' @export rf_est_FS
rf_est_FS <- function(init_r=0.01, epsilon=0.001, depth_Ref, depth_Alt, OPGP,
                      sexSpec=F, trace=F, noFam=1, method = "EM", ...){
  
  ## Do some checks
  if(!is.list(depth_Ref) | !is.list(depth_Alt) | !is.list(OPGP))
    stop("Arguments for read count matrices and vector of OPGPs are required to be list objects")
  if( !is.numeric(noFam) || noFam < 1 || noFam != round(noFam) || !is.finite(noFam))
    stop("The number of families needs to be a finite positive number")
  if(noFam != length(depth_Ref) | noFam != length(depth_Alt) | noFam != length(OPGP) )
    stop("The number of read count matrices or OPGP vectors do not match the number of families specified")
  if( !is.null(init_r) & !is.numeric(init_r) )
    stop("Starting values for the recombination fraction needs to be a numeric vector or integer or a NULL object")
  if( (!is.null(epsilon) & !is.numeric(epsilon)) || (!is.null(epsilon) && (epsilon <= 0 | epsilon >= 1)) )
    stop("Starting values for the error parameters needs to be a single numeric value in the interval (0,1) or a NULL object")
  if( !is.logical(trace) || is.na(trace) )
    trace = FALSE
  if( !is.logical(sexSpec) || is.na(sexSpec) )
    sexSpec = FALSE
  if(!(method %in% c("EM","optim")))
     stop("Specified optimization method is unknown. Please select one of 'EM' or 'optim'")

  ## Check the read count matrices
  if(any(unlist(lapply(depth_Ref,function(x) !is.numeric(x) || any( x<0 | !is.finite(x)) || any(!(x == round(x)))))))
    stop("At least one read count matrix for the reference allele is missing or invalid")
  if(any(unlist(lapply(depth_Alt,function(x) !is.numeric(x) || any( x<0 | !is.finite(x)) || any(!(x == round(x)))))))
    stop("At least one read count matrix for the alternate allele is missing or invalid")
  if(any(unlist(lapply(OPGP, function(x) !is.numeric(x) || !is.vector(x) || any(!(x %in% 1:16)) ))))
    stop("At least OPGP vector is missing or invalid")
     
  nInd <- lapply(depth_Ref,nrow)  # number of individuals
  nSnps <- ncol(depth_Ref[[1]])   # number of SNPs
  
  ## check inputs are of required type for C functions
  if(!is.numeric(init_r)|is.integer(init_r))
    init_r <- as.numeric(init_r)
  for(fam in 1:noFam){
    if(!is.integer(depth_Ref[[fam]]))
      depth_Ref[[fam]] <- matrix(as.integer(depth_Ref[[fam]]), nrow=nInd[[fam]], ncol=nSnps)
    if(!is.integer(depth_Alt[[fam]]))
      depth_Alt[[fam]] <- matrix(as.integer(depth_Alt[[fam]]), nrow=nInd[[fam]], ncol=nSnps)
    if(!is.integer(OPGP[[fam]]))
      OPGP[[fam]] <- as.integer(OPGP[[fam]])
  }
  if(!is.integer(noFam))
    noFam <- as.integer(noFam)
  
  if(method=="optim"){
  
    # Arguments for the optim function
    optim.arg <- list(...)
    if(length(optim.arg) == 0)
      optim.arg <- list(maxit = 1000, reltol=1e-15)
    
    ## Compute the K matrix for heterozygous genotypes
    bcoef_mat <- Kab <- vector(mode="list", length=noFam)
    for(fam in 1:noFam){
      bcoef_mat[[fam]] <- choose(depth_Ref[[fam]]+depth_Alt[[fam]],depth_Ref[[fam]])
      Kab[[fam]] <- bcoef_mat[[fam]]*(1/2)^(depth_Ref[[fam]]+depth_Alt[[fam]])
    }
    
    ## If we want to estimate sex-specific r.f.'s
    if(sexSpec){
      
      # Work out the indices of the r.f. parameters of each sex
      ps <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% 1:8)))))[-1] - 1
      ms <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% c(1:4,9:12))))))[-1] - 1
      npar <- c(length(ps),length(ms))
      
      # Determine the initial values
      if(length(init_r)==1) 
        para <- logit2(rep(init_r,sum(npar)))
      else if(length(init_r) != sum(npar)) 
        para <- logit2(rep(0.1,sum(npar)))
      else
        para <- init_r
      # sequencing error
      if(length(epsilon) != 1 & !is.null(epsilon))
        para <- c(para,logit(0.001))
      else if(!is.null(epsilon))
        para <- c(para,logit(epsilon))
    
      ## Are we estimating the error parameters?
      seqErr=!is.null(epsilon)
      
      ## Find MLE
      optim.MLE <- optim(para,ll_fs_ss_mp_scaled_err,method="BFGS",control=optim.arg,
                         depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
                         nInd=nInd,nSnps=nSnps,OPGP=OPGP,ps=ps,ms=ms,npar=npar,noFam=noFam,
                         seqErr=!is.null(epsilon))
    }
    else{
      # Determine the initial values
      if(length(init_r)==1) 
        para <- logit2(rep(init_r,nSnps-1))
      else if(length(init_r) != nSnps-1) 
        para <- logit2(rep(0.1,nSnps-1))
      else
        para <- init_r
      # sequencing error
      if(length(epsilon) != 1 & !is.null(epsilon))
        para <- c(para,logit(0.001))
      else if(!is.null(epsilon))
        para <- c(para,logit(epsilon))
      
      ## Are we estimating the error parameters?
      seqErr=!is.null(epsilon)
      
      ## Find MLE
      optim.MLE <- optim(para,ll_fs_mp_scaled_err,method="BFGS",control=optim.arg,
                         depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
                         nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam,
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
    if(sexSpec)
      return(list(rf_p=inv.logit2(optim.MLE$par[1:npar[1]]),rf_m=inv.logit2(optim.MLE$par[npar[1]+1:npar[2]]),
                  epsilon=ifelse(seqErr,inv.logit(optim.MLE$par[sum(npar)+1]),0),
                  loglik=optim.MLE$value))
    else
      return(list(rf=inv.logit2(optim.MLE$par[1:(nSnps-1)]), 
                  epsilon=ifelse(seqErr,inv.logit(optim.MLE$par[nSnps]),0),
                  loglik=optim.MLE$value))
  }
  else{ # EM algorithm approach
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

    # sequencing error
    if(length(epsilon) != 1 & !is.null(epsilon))
      epsilon <- 0.001
    
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
    depth_Ref_mat = do.call(what = "rbind",depth_Ref)
    depth_Alt_mat = do.call(what = "rbind",depth_Alt)
    
    ## Are we estimating the error parameters?
    seqErr=!is.null(epsilon)
    if(is.null(epsilon))
      epsilon = 0
    
    EMout <- .Call("EM_HMM", init_r, epsilon, depth_Ref_mat, depth_Alt_mat, OPGPmat,
                   noFam, unlist(nInd), nSnps, sexSpec, seqErr, EM.arg, as.integer(ss_rf))

    if(sexSpec){
      return(list(rf_p=EMout[[1]][ps],rf_m=EMout[[1]][nSnps-1+ms],
                  epsilon=EMout[[2]],
                  loglik=EMout[[3]]))
    }
    else
      return(list(rf=EMout[[1]][1:(nSnps-1)], 
                  epsilon=EMout[[2]],
                  loglik=EMout[[3]]))
    
  }
}


## recombination estimates for case where the phase is unkonwn.
## The r.f.'s are sex-specific and constrained to the range [0,1]
rf_est_FS_UP <- function(depth_Ref, depth_Alt, config, epsilon, method="EM", trace=F, ...){
  
  ## Check imputs
  if( any( depth_Ref<0 | !is.finite(depth_Ref)) || any(!(depth_Ref == round(depth_Ref))))
    stop("The read count matrix for the reference allele is invalid")
  if( any( depth_Alt<0 | !is.finite(depth_Alt)) || any(!(depth_Alt == round(depth_Alt))))
    stop("The read count matrix for the alternate allele is invalid")
  if( !is.vector(config) || !is.numeric(config) || any(!(config %in% 1:9)) )
    stop("Invalid config vector. It must be a numeric vector with entries between 1 and 9.")
  if( !is.logical(trace) || is.na(trace) )
    trace = FALSE
  if(!(method %in% c("EM","optim")))
     stop("Specified optimization method is unknown. Please select one of 'EM' or 'optim'")
  
  nInd <- nrow(depth_Ref)  # number of individuals
  nSnps <- ncol(depth_Ref)  # number of SNPs
  if(!is.integer(depth_Ref))
    depth_Ref <- matrix(as.integer(depth_Ref), nrow=nInd, ncol=nSnps)
  if(!is.integer(depth_Alt))
     depth_Alt <- matrix(as.integer(depth_Alt), nrow=nInd, ncol=nSnps)
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
    bcoef_mat <- choose(depth_Ref+depth_Alt,depth_Ref)
    Kab <- bcoef_mat*(1/2)^(depth_Ref+depth_Alt)
    
    ## Are we estimating the error parameters?
    seqErr=!is.null(epsilon)
    
    para <- logit(rep(0.5,sum(npar)))
    # sequencing error
    if(length(epsilon) != 1 & !is.null(epsilon))
      para <- c(para,logit(0.01))
    else if(!is.null(epsilon))
      para <- c(para,logit(epsilon))
    
    if(nSnps > 2){
      ## Find MLE
      optim.MLE <- optim(para,ll_fs_up_ss_scaled_err,method="BFGS",control=optim.arg,
                           depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
                         nInd=nInd,nSnps=nSnps,config=config,ps=ps,ms=ms,npar=npar,
                         seqErr=seqErr)
      # Print out the output from the optim procedure (if specified)
      if(trace){
        print(optim.MLE)
      }
      
      # Check for convergence
      if(optim.MLE$convergence != 0)
        warning(paste0('Optimization failed to converge properly with error ',optim.MLE$convergence,'\n smallest MLE estimate is: ', round(min(optim.MLE$par),6)))
      
      # Return the MLEs
      return(list(rf_p=inv.logit(optim.MLE$par[1:npar[1]]),rf_m=inv.logit(optim.MLE$par[npar[1]+1:npar[2]]),
                  epsilon=inv.logit(optim.MLE$par[sum(npar)+1])))
    } 
    else if(nSnps == 2){
      ## If both SNPs are informative, need to use the Nelder-Mead to distinguish between the two sexes.
      if(all(config == 1)){
        optim.MLE <- optim(para,ll_fs_up_ss_scaled_err,method="Nelder-Mead",control=optim.arg,
                           depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
                           nInd=nInd,nSnps=nSnps,config=config,ps=ps,ms=ms,npar=npar,
                           seqErr=seqErr)
      }
      ## Otherwise, proceed as normal
      else{
        optim.MLE <- optim(para,ll_fs_up_ss_scaled_err,method="BFGS",control=optim.arg,
                           depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
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
                  epsilon=inv.logit(optim.MLE$par[sum(npar)+1])))
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
    seqErr=!is.null(epsilon)
    
    EMout <- .Call("EM_HMM_UP", rep(0.5,(nSnps-1)*2), epsilon, depth_Ref, depth_Alt, config,
          as.integer(1), nInd, nSnps, seqErr, EM.arg, as.integer(ss_rf))
    return(list(rf_p=EMout[[1]][ps],rf_m=EMout[[1]][nSnps-1+ms],
                epsilon=EMout[[2]],
                loglik=EMout[[3]]))
  }
}



