### R Script of functions for computing the recombination fraction estimate 
### for linkage analysis using GBS data.
### Author: Timothy Bilton
### Date: 18/01/17
### Edited: 03/04/17


## Function for computing the recombination fraction when the parental phase is known
rf_est_FS <- function(init_r=NULL, genon, depth, OPGP, sexSpec=F, trace=F, noFam=1, ...){
  
  if((is.list(genon) & noFam != length(genon))|
     (is.list(depth) & noFam != length(depth))|
     (is.list(OPGP) & noFam != length(OPGP)))
    stop("Error: The number of genon and depth matrices given do not match the number of families expected")
  
  else{
    # Arguments for the optim function
    optim.arg <- list(...)
    if(length(optim.arg) == 0)
      optim.arg <- list(maxit = 1000, reltol=1e-10)
    
    ## if only one family.
    if(noFam == 1){
    
      ## Convert the geno data into alternative form
      genon <- abs(2-genon+1)
      genon[which(is.na(genon))] <- 4
    
      nInd <- nrow(genon)   # number of individuals
      nSnps <- ncol(genon)  # number of SNPs
      
      ## check inputs are of required type for C functions
      if(!is.numeric(init_r)|is.integer(init_r))
        init_r <- as.numeric(init_r)
      if(!is.numeric(genon)|is.integer(genon))
        genon <- genon + 0
      if(!is.numeric(depth)|is.integer(depth))
        depth <- depth + 0
      if(!is.numeric(OPGP)|is.integer(OPGP))
        OPGP <- as.numeric(OPGP)
    
      if(sexSpec){
    
        # Work out the indices of the r.f. parameters of each sex
        ps <- which(OPGP %in% 1:6)[-1] - 1
        ms <- which(OPGP %in% c(1:4,7:8))[-1] - 1
        npar <- c(length(ps),length(ms))
        
        # Determine the initial values
        if((length(init_r)==1) & (is.numeric(init_r))) init_r <- rep(init_r,sum(npar))
        else if((length(init_r) != sum(npar)) | !(is.numeric(init_r))) init_r <- rep(0.1,sum(npar))
        
        ## Find MLE
        if(scaled){
          optim.MLE <- optim(logit2(init_r),ll_fs_ss_scaled,method="BFGS",control=optim.arg,
                             genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,OPGP=OPGP,ps=ps,ms=ms,npar=npar, hessian=hessian)
        }
        else{
          optim.MLE <- optim(logit2(init_r),ll_fs_ss,method="BFGS",control=optim.arg,
                             genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,OPGP=OPGP,ps=ps,ms=ms,npar=npar, hessian=hessian)
        }
      }
      else{
        # Determine the initial values
        if((length(init_r)==1) & (is.numeric(init_r))) init_r <- rep(init_r,nSnps-1)
        else if((length(init_r) != nSnps) & (!is.numeric(init_r))) init_r <- rep(0.1,nSnps-1)
        
        ## Find MLE
        if(scaled){
          optim.MLE <- optim(logit2(init_r),ll_fs_scaled,method="BFGS",control=optim.arg,
                           genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,OPGP=OPGP, hessian=hessian)
        } else{
          optim.MLE <- optim(logit2(init_r),ll_fs,method="BFGS",control=optim.arg,
                             genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,OPGP=OPGP, hessian=hessian)
        }
      }
    }
    else{
      ## Convert the geno data into alternative form
      genon <- lapply(genon,function(x){
        x <- abs(2-x+1)
        x[which(is.na(x))] <- 4
        return(x)
      })
      
      nInd <- lapply(genon,nrow)  # number of individuals
      nSnps <- ncol(genon[[1]]) # number of SNPs
      
      ## check inputs are of required type for C functions
      if(!is.numeric(init_r)|is.integer(init_r))
        init_r <- as.numeric(init_r)
      for(fam in 1:noFam){
        if(!is.numeric(genon[[fam]])|is.integer(genon[[fam]]))
          genon[[fam]] <- genon[[fam]] + 0
        if(!is.numeric(depth[[fam]])|is.integer(depth[[fam]]))
          depth[[fam]] <- depth[[fam]] + 0
        if(!is.numeric(OPGP[[fam]])|is.integer(OPGP[[fam]]))
          OPGP <- as.numeric(OPGP[[fam]])
      }
      
      if(sexSpec){
        
        # Work out the indices of the r.f. parameters of each sex
        ps <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% 1:6)))))[-1] - 1
        ms <- sort(unique(unlist(lapply(OPGP,function(x) which(x %in% c(1:4,7:8))))))[-1] - 1
        npar <- c(length(ps),length(ms))
        
        # Determine the initial values
        if((length(init_r)==1) & (is.numeric(init_r))) init_r <- rep(init_r,sum(npar))
        else if((length(init_r) != sum(npar)) | !(is.numeric(init_r))) init_r <- rep(0.1,sum(npar))
        
        ## Find MLE
        if(scaled){
          optim.MLE <- optim(logit2(init_r),ll_fs_ss_mp_scaled,method="BFGS",control=optim.arg,
                             genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,OPGP=OPGP,ps=ps,ms=ms,npar=npar, noFam=noFam, hessian=hessian)
        }
        else{
          optim.MLE <- optim(logit2(init_r),ll_fs_ss_mp,method="BFGS",control=optim.arg,
                             genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,OPGP=OPGP,ps=ps,ms=ms,npar=npar, noFam=noFam, hessian=hessian)
        }
      }
      else{
        # Determine the initial values
        if((length(init_r)==1) & (is.numeric(init_r))) init_r <- rep(init_r,nSnps-1)
        else if((length(init_r) != nSnps) & (!is.numeric(init_r))) init_r <- rep(0.1,nSnps-1)
        
        ## Find MLE
        if(scaled){
          optim.MLE <- optim(logit2(init_r),ll_fs_mp_scaled,method="BFGS",control=optim.arg,
                             genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,OPGP=OPGP, noFam=noFam, hessian=hessian)
        } else{
          optim.MLE <- optim(logit2(init_r),ll_fs_mp,method="BFGS",control=optim.arg,
                             genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,OPGP=OPGP, noFam=noFam, hessian=hessian)
        }
      }
    }
    # Print out the output from the optim procedure (if specified)
    if(trace){
      print(optim.MLE)
    }
    # Check for convergence
    if(optim.MLE$convergence != 0)
      warning(paste0('Optimization failed to converge properly with error ',optim.MLE$convergence,'\n smallest MLE estimate is: ', round(min(optim.MLE$par),6)))
    # Return the MLEs
    if(sexSpec){
      if(hessian)
        return(list(inv.logit2(optim.MLE$par[1:npar[1]]),inv.logit2(optim.MLE$par[npar[1]+1:npar[2]]), optim.MLE$value,optim.MLE$hessian))
      else
        return(list(inv.logit2(optim.MLE$par[1:npar[1]]),inv.logit2(optim.MLE$par[npar[1]+1:npar[2]]), optim.MLE$value))
    }
    else{
      if(hessian)
        return(list(inv.logit2(optim.MLE$par), optim.MLE$value, optim.MLE$hessian))
      else
        return(list(inv.logit2(optim.MLE$par), optim.MLE$value))
    }
  }
}


#### Group 2: Functions required for finding the MLE of the GBS HMM where OPGP is not known #####
## The r.f. is now on the range [0,1] and we assume the default phase is coupling
## so that r > 0.5 impies the repulsion phase. 

#### Function for finding the MLE of the GBS HMM
# Input Variables
#  - ini_r_f: Initial values for the vector of the paternal r.f.'s 
#  - ini_r_m: Initial values for the vector of the maternal r.f.'s 
#  - genon: The matrix of genotype calls. 0=BB, 1=AB, 2=AA, NA=missing
#  - depth: The matrix of depth calls. Positive for non-missing genotype calls.
#  - config: Segregation type of each marker. 1=informative, 2=paternal segregating, 3=maternal segregating
#  - trace: If TRUE prints out the output for optim
rf_est_FS_UP <- function(genon, depth, config, trace=F, scaled=T, ...){
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
  if(!is.numeric(genon)|is.integer(genon))
    genon <- genon + 0
  if(!is.numeric(depth)|is.integer(depth))
    depth <- depth + 0
  if(!is.numeric(config)|is.integer(config))
    config <- as.numeric(config)
  
  # Work out the indices of the r.f. parameters of each sex
  ps <- which(config %in% c(1,2))[-1] - 1
  ms <- which(config %in% c(1,3))[-1] - 1
  npar <- c(length(ps),length(ms))
  
  if(nSnps > 2){
    ## Find MLE
    if(scaled){
      optim.MLE <- optim(logit(rep(0.5,sum(npar))),ll_fs_up_ss_scaled,method="BFGS",control=optim.arg,
                         genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,config=config,ps=ps,ms=ms,npar=npar)
    } else{
      optim.MLE <- optim(logit(rep(0.5,sum(npar))),ll_fs_up_ss,method="BFGS",control=optim.arg,
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
    return(list(inv.logit(optim.MLE$par[1:npar[1]]),inv.logit(optim.MLE$par[npar[1]+1:npar[2]])))
  } 
  else if(nSnps == 2){
    ## If both SNPs are informative, need to use the Nelder-Mead to distinguish between the two sexes.
    if(all(config == 1)){
      optim.MLE <- optim(logit(rep(0.5,sum(npar))),ll_fs_up_ss,method="Nelder-Mead",control=optim.arg,
                         genon=genon,depth=depth,nInd=nInd,nSnps=nSnps,config=config,ps=ps,ms=ms,npar=npar)
    }
    ## Otherwise, proceed as normal
    else{
      optim.MLE <- optim(logit(rep(0.5,sum(npar))),ll_fs_up_ss,method="BFGS",control=optim.arg,
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
    return(list(inv.logit(optim.MLE$par[npar[1]]),inv.logit(optim.MLE$par[npar[1]+npar[2]])))
  }
}



