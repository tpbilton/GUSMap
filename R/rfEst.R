### R Script of functions for computing the recombination fraction estimate 
### for linkage analysis using sequencing data.
### Author: Timothy Bilton
### Date: 18/01/17
### Edited: 19/06/17


## Function for computing the recombination fraction when the parental phase is known
rf_est_FS <- function(init_r=NULL, genon, depth, OPGP,
                      sexSpec=F, trace=F, noFam=1, ...){
  
  ## Do some checks
  if(!is.list(genon) | !is.list(depth) | !is.list(OPGP))
    stop("Arguments for genon, depth and OPGP are required to be lists objects")
  if( !is.numeric(noFam) || noFam < 1 || noFam != round(noFam) || !is.finite(noFam))
    stop("The number of families needs to be a finite positive number")
  if(noFam != length(genon) | noFam != length(depth) | noFam != length(OPGP) )
    stop("The number of genon and depth matrices or OPGP vectors do not match the number of families specified")
  if( !is.null(init_r) & !is.numeric(init_r) )
    stop("Starting values for the recombination fraction needs to be a numeric vector or integer or a NULL object")
  if( !is.logical(trace) || is.na(trace) )
    trace = FALSE
  if( !is.logical(sexSpec) || is.na(sexSpec) )
    sexSpec = FALSE

  # Arguments for the optim function
  optim.arg <- list(...)
  if(length(optim.arg) == 0)
    optim.arg <- list(maxit = 1000, reltol=1e-10)
  
  ## Convert the geno data into alternative form
  genon <- lapply(genon,function(x){
    if(!is.matrix(x) || !is.numeric(x) || any(!(x %in%  0:2 | is.na(x))))
      return(NULL)
    else{
      x <- abs(2-x+1)
      x[which(is.na(x))] <- 4
      return(x)
    }
  })
  ## Check the genon and depth matrices
  if(any(unlist(lapply(genon,is.null))))
    stop("At least one genon matrix is missing or invalid")
  if(any(unlist(lapply(depth,function(x) !is.numeric(x) || any( x<0 | !is.finite(x)) || any(!(x == round(x)))))))
    stop("At least one depth matrix is missing or invalid")
  if(any(unlist(lapply(OPGP, function(x) !is.numeric(x) || !is.vector(x) || any(!(x %in% 1:9)) ))))
    stop("At least OPGP vector is missing or invalid")
     
  nInd <- lapply(genon,nrow)  # number of individuals
  nSnps <- ncol(genon[[1]])   # number of SNPs
  
  ## check inputs are of required type for C functions
  if(!is.numeric(init_r)|is.integer(init_r))
    init_r <- as.numeric(init_r)
  for(fam in 1:noFam){
    if(is.integer(genon[[fam]]))
      genon[[fam]] <- genon[[fam]] + 0
    if(is.integer(depth[[fam]]))
      depth[[fam]] <- depth[[fam]] + 0
    if(is.integer(OPGP[[fam]]))
      OPGP <- as.numeric(OPGP[[fam]])
  }
  
  ## If we want to estimate sex-specific r.f.'s
  if(sexSpec){
    
    # Work out the indices of the r.f. parameters of each sex
    ps <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% 1:6)))))[-1] - 1
    ms <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% c(1:4,7:8))))))[-1] - 1
    npar <- c(length(ps),length(ms))
    
    # Determine the initial values
    if((length(init_r)==1) & (is.numeric(init_r))) init_r <- rep(init_r,sum(npar))
    else if((length(init_r) != sum(npar)) | !(is.numeric(init_r))) init_r <- rep(0.1,sum(npar))
    
    ## Find MLE
    optim.MLE <- optim(logit2(init_r),ll_fs_ss_mp_scaled,method="BFGS",control=optim.arg,
                         genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,OPGP=OPGP,ps=ps,ms=ms,npar=npar,noFam=noFam)
  }
  else{
    # Determine the initial values
    if((length(init_r)==1) & (is.numeric(init_r))) init_r <- rep(init_r,nSnps-1)
    else if((length(init_r) != nSnps) | (!is.numeric(init_r))) init_r <- rep(0.1,nSnps-1)
    
    ## Find MLE
    optim.MLE <- optim(logit2(init_r),ll_fs_mp_scaled,method="BFGS",control=optim.arg,
                         genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,OPGP=OPGP,noFam=noFam)
  }

  # Print out the output from the optim procedure (if specified)
  if(trace){
    print(optim.MLE)
  }
  # Check for convergence
  if(optim.MLE$convergence != 0)
    warning(paste0('Optimization failed to converge properly with error ',optim.MLE$convergence,'\n smallest MLE estimate is: ', round(min(optim.MLE$par),6)))
  # Return the MLEs
  if(sexSpec)
    return(list(rf_p=inv.logit2(optim.MLE$par[1:npar[1]]),rf_m=inv.logit2(optim.MLE$par[npar[1]+1:npar[2]]), loglik=optim.MLE$value))
  else
    return(list(rf=inv.logit2(optim.MLE$par), loglik=optim.MLE$value))
}


## recombination estimates for case where the phase is unkonwn.
## The r.f.'s are sex-specific and constrained to the range [0,1]
rf_est_FS_UP <- function(genon, depth, config, trace=F, ...){
  
  ## Check imputs
  if( !is.matrix(genon) || !is.numeric(genon) || any(!(genon %in% 0:2 | is.na(genon))) )
    stop("Invalid genon matrix. It must be numeric and have entries of 0, 1, 2 or NA")
  if( !is.matrix(depth) || !is.numeric(depth) || any(depth<0 || !is.finite(depth)) || any(!(depth == round(depth))) )
    stop("Invalid depth matrix. It must be numeric with entries that are finite positive integers")
  if( !is.vector(config) || !is.numeric(config) || any(!(config %in% 1:4)) )
    stop("Invalid config vector. It must be a numeric vector with entries 1, 2, 3 or 4.")
  if( !is.logical(trace) || is.na(trace) )
    trace = FALSE
  
  ## Convert the geno data into alternative form to allow indexing of lists
  genon <- abs(2-genon+1)
  genon[which(is.na(genon))] <- 4
  
  nInd <- nrow(genon)  # number of individuals
  nSnps <- ncol(genon)  # number of SNPs
  
  # Arguments for the optim function
  optim.arg <- list(...)
  if(length(optim.arg) == 0)
    optim.arg <- list(maxit = 1000, reltol=1e-10)
  
  ## check inputs are of required type for C functions
  if(is.integer(genon))
    genon <- genon + 0
  if(is.integer(depth))
    depth <- depth + 0
  if(is.integer(config))
    config <- as.numeric(config)
  
  # Work out the indices of the r.f. parameters of each sex
  ps <- which(config %in% c(1,2))[-1] - 1
  ms <- which(config %in% c(1,3))[-1] - 1
  npar <- c(length(ps),length(ms))
  
  if(nSnps > 2){
    ## Find MLE
    optim.MLE <- optim(logit(rep(0.5,sum(npar))),ll_fs_up_ss_scaled,method="BFGS",control=optim.arg,
                         genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,config=config,ps=ps,ms=ms,npar=npar)
    # Print out the output from the optim procedure (if specified)
    if(trace){
      print(optim.MLE)
    }
    
    # Check for convergence
    if(optim.MLE$convergence != 0)
      warning(paste0('Optimization failed to converge properly with error ',optim.MLE$convergence,'\n smallest MLE estimate is: ', round(min(optim.MLE$par),6)))
    
    # Return the MLEs
    return(list(rf_p=inv.logit(optim.MLE$par[1:npar[1]]),rf_m=inv.logit(optim.MLE$par[npar[1]+1:npar[2]])))
  } 
  else if(nSnps == 2){
    ## If both SNPs are informative, need to use the Nelder-Mead to distinguish between the two sexes.
    if(all(config == 1)){
      optim.MLE <- optim(logit(rep(0.5,sum(npar))),ll_fs_up_ss_scaled,method="Nelder-Mead",control=optim.arg,
                         genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,config=config,ps=ps,ms=ms,npar=npar)
    }
    ## Otherwise, proceed as normal
    else{
      optim.MLE <- optim(logit(rep(0.5,sum(npar))),ll_fs_up_ss_scaled,method="BFGS",control=optim.arg,
                         genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,config=config,ps=ps,ms=ms,npar=npar)
    }
    # Print out the output from the optim procedure (if specified)
    if(trace){
      print(optim.MLE)
    }
    
    # Check for convergence
    if(optim.MLE$convergence != 0)
      warning(paste0('Optimization failed to converge properly with error ',optim.MLE$convergence,'\n smallest MLE estimate is: ', round(min(optim.MLE$par),6)))
    
    # Return the MLEs
    return(list(rf_p=inv.logit(optim.MLE$par[npar[1]]),rf_m=inv.logit(optim.MLE$par[npar[1]+npar[2]])))
  }
}



