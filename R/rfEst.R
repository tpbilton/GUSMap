### R Script of functions for computing the recombination fraction estimate 
### for linkage analysis using sequencing data.
### Author: Timothy Bilton
### Date: 18/01/17
### Edited: 19/06/17


## Function for computing the recombination fraction when the parental phase is known
rf_est_FS <- function(init_r=0.01, epsilon=0.001, depth_Ref, depth_Alt, OPGP,
                      sexSpec=F, trace=F, noFam=1, ...){
  
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

  # Arguments for the optim function
  optim.arg <- list(...)
  if(length(optim.arg) == 0)
    optim.arg <- list(maxit = 1000, reltol=1e-15)
  
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
    if(is.integer(OPGP[[fam]]))
      OPGP <- as.numeric(OPGP[[fam]])
  }
  
  ## Compute the K matrix for heterozygous genotypes
  bcoef_mat <- Kab <- vector(mode="list", length=noFam)
  for(fam in 1:noFam){
    bcoef_mat[[fam]] <- choose(depth_Ref[[fam]]+depth_Alt[[fam]],depth_Ref[[fam]])
    Kab[[fam]] <- bcoef_mat[[fam]]*(1/2)^(depth_Ref[[fam]]+depth_Alt[[fam]])
  }
  
  ## If multi-point
  if(nSnps > 2){
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
  else if (nSnps==2){
    
    para <- logit2(0.2)
    # sequencing error
    if(length(epsilon) != 1 & !is.null(epsilon))
      para <- c(para,logit(0.001))
    else if(!is.null(epsilon))
      para <- c(para,logit(epsilon))
    
    ## Are we estimating the error parameters?
    seqErr=!is.null(epsilon)
    
    ## Find MLE
    optim.MLE <- optim(para,ll_fs_mp_scaled_err,method="BFGS",control=list(reltol=1e-10),
                       depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
                       nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam,
                       seqErr=seqErr)
    
    ## Compute the log score
    LOD <-  -(optim.MLE$value - ll_fs_mp_scaled_err( ifelse(seqErr,c(1000,optim.MLE$par[2]), 1000),
                                                   depth_Ref=depth_Ref,depth_Alt=depth_Alt,bcoef_mat=bcoef_mat,Kab=Kab,
                                                   nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam,
                                                   seqErr=seqErr))
    
    # Print out the output from the optim procedure (if specified)
    if(trace){
      print(optim.MLE)
    }
    
    return(list(rf=inv.logit2(optim.MLE$par[1]), 
                epsilon=0,
                loglik=optim.MLE$value,LOD=LOD))
  }
  
}


## recombination estimates for case where the phase is unkonwn.
## The r.f.'s are sex-specific and constrained to the range [0,1]
rf_est_FS_UP <- function(depth_Ref, depth_Alt, config, epsilon, trace=F, ...){
  
  ## Check imputs
  if( any( depth_Ref<0 | !is.finite(depth_Ref)) || any(!(depth_Ref == round(depth_Ref))))
    stop("The read count matrix for the reference allele is invalid")
  if( any( depth_Alt<0 | !is.finite(depth_Alt)) || any(!(depth_Alt == round(depth_Alt))))
    stop("The read count matrix for the alternate allele is invalid")
  if( !is.vector(config) || !is.numeric(config) || any(!(config %in% 1:9)) )
    stop("Invalid config vector. It must be a numeric vector with entries between 1 and 9.")
  if( !is.logical(trace) || is.na(trace) )
    trace = FALSE
  
  nInd <- nrow(depth_Ref)  # number of individuals
  nSnps <- ncol(depth_Ref)  # number of SNPs
  
  # Arguments for the optim function
  optim.arg <- list(...)
  if(length(optim.arg) == 0)
    optim.arg <- list(maxit = 1000, reltol=1e-10)
  
  ## check inputs are of required type for C functions
  if(is.integer(config))
    config <- as.numeric(config)
  
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



