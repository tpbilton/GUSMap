

#### Function for creating a particular population structure
createPop <- function(R6obj, pop = c("full-sib"), ...){
  
  ## Make new R6 object depending on the family type.
  if(pop == "full-sib")
    newObj <- FS$new(R6obj)
  else
    stop(paste("Population structure",pop,"has not yet be implemented\n"))
  
  ## Create the population
  return(makePop(newObj,...))
}



### Generic method for creating a population
makePop <- function(obj, ...){
  UseMethod("makePop")
}

### Make a full-sib family population
makePop.FS <- function(R6obj, famInfo, filter=list(MAF=0.05, MISS=0.2, BIN=0, DEPTH=6, PVALUE=0.05)){
  
  ## Do some checks
  if(is.null(filter$MAF) || filter$MAF<0 || filter$MAF>1 || !is.numeric(filter$MAF))
    stop("Minor allele frequency filter has not be specifies or is invalid.")
  
  ## Define variables that will be used.
  noFam <- length(famInfo)
  config_all <- config_infer_all <- nSnps_all <- nInd_all <- indx <- indID_all <- vector(mode = "list", length = noFam)
  
  cat("-------------\n")
  cat("Processing Data.\n\n")
  
  cat("Filtering criteria for removing SNPs :\n")
  cat("Minor allele frequency (MAF) < ", filter$MAF,"\n")
  cat("Percentage of missing genotypes > ", filter$MISS*100,"%\n\n",sep="")
  
  ## Extract the private variables we want
  indID <- R6obj$.__enclos_env__$private$indID
  nSnps <- R6obj$.__enclos_env__$private$nSnps
  
  ## extract the data and format correct for each family.
  for(fam in 1:noFam){
    cat("Processing Family ",names(famInfo)[fam],".\n\n",sep="")
    mum <- famInfo[[fam]]$parents$Mother
    dad <- famInfo[[fam]]$parents$Father
    ## index the parents
    mumIndx <- which(indID %in% mum)
    if(length(mumIndx) == 0)
      stop(paste0("Mother ID not found family ",fam,"."))
    dadIndx <- which(indID %in% dad)
    if(length(dadIndx) == 0)
      stop(paste0("Father ID not found family ",fam,"."))
    ## index the progeny
    progIndx <- which(indID %in% famInfo[[fam]]$progeny)
    nInd <- length(progIndx)
    indID_all[[fam]] <- indID[progIndx]
    ## Subset the genon and depth matrices
    genon     <- R6obj$.__enclos_env__$private$genon[progIndx,]
    depth_Ref <- R6obj$.__enclos_env__$private$depth_Ref[progIndx,]
    depth_Alt <- R6obj$.__enclos_env__$private$depth_Alt[progIndx,]

    ## Determine the segregation types of the loci
    genon_mum <- matrix(R6obj$.__enclos_env__$private$genon[mumIndx,], nrow=length(mumIndx), ncol=nSnps) 
    genon_dad <- matrix(R6obj$.__enclos_env__$private$genon[dadIndx,], nrow=length(mumIndx), ncol=nSnps)
    depth_mum <- matrix(R6obj$.__enclos_env__$private$depth_Ref[mumIndx,] +
                          R6obj$.__enclos_env__$private$depth_Alt[mumIndx,], nrow=length(mumIndx), ncol=nSnps)
    depth_dad <- matrix(R6obj$.__enclos_env__$private$depth_Ref[dadIndx,] +
                          R6obj$.__enclos_env__$private$depth_Alt[dadIndx,], nrow=length(mumIndx), ncol=nSnps)
    
    ## Determine segregation type of each SNP if possible
    config <- unlist(sapply(1:nSnps,function(x){
      x_p = genon_dad[,x]; x_m = genon_mum[,x]
      d_p = depth_dad[,x]; d_m = depth_mum[,x]
      if(sum(d_p) > filter$DEPTH & sum(d_p) > filter$DEPTH ){
        if(any(x_p==1,na.rm=T) & any(x_m==1,na.rm=T))
          return(1)
        else if(all(x_m==2,na.rm=T) & (any(x_p==1, na.rm=T) | all(x_p %in% c(0,2), na.rm=T)))
          return(2)
        else if(all(x_m==0,na.rm=T) & (any(x_p==1, na.rm=T) | all(x_p %in% c(0,2), na.rm=T)))
          return(3)
        else if(all(x_p==2,na.rm=T) & (any(x_m==1, na.rm=T) | all(x_m %in% c(0,2), na.rm=T)))
          return(4)
        else if(all(x_p==0,na.rm=T) & (any(x_m==1, na.rm=T) | all(x_m %in% c(0,2), na.rm=T)))
          return(5)
        else
          return(NA)
      }
      else return(NA)
    }))
    
    #### Segregation test to determine if the SNPs have been miss-classified
    seg_Dis <- sapply(1:nSnps,function(x){
      if(is.na(config[x]))
        return(NA)
      else{
        d = depth_Ref[,x] + depth_Alt[,x]
        g = genon[,x]
        K = sum(1/2^(d[which(d != 0)])*0.5)/sum(d != 0)
        nAA = sum(g==2, na.rm=T)
        nAB = sum(g==1, na.rm=T)
        nBB = sum(g==0, na.rm=T)
        ## check that there are sufficient data to perform the chisq test
        if(sum(nAA+nAB+nBB)/length(g) <= (1-filter$MISS))
          return(NA)
        else if(config[x] == 1){
          exp_prob <- c(0.25 + K,0.5 - 2*K, 0.25 + K)
          ctest <- chisq.test(c(nBB,nAB,nAA), p = exp_prob)
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
        else if(config[x] %in% c(2,4)){
          exp_prob <- c(K, 0.5 - 2*K, 0.5 + K)
          ctest <- chisq.test(c(nBB,nAB,nAA), p = exp_prob)
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
        else if(config[x] %in% c(3,5)){
          exp_prob <- c(0.5 + K, 0.5 - 2*K, K)
          ctest <- chisq.test(c(nBB,nAB,nAA), p = exp_prob)
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
      }
    },simplify = T)
    config[which(seg_Dis)] <- NA
    
    ## Run the filtering of the progeny SNPs
    MAF <- colMeans(genon, na.rm=T)/2
    MAF <- pmin(MAF,1-MAF)
    miss <- apply(genon,2, function(x) sum(is.na(x))/length(x))
    
    ## Infer geotypes for over SNPs that have passed the MAF and MISS thresholds
    #propHeter <- sapply(1:nSnps, function(x) sum(genon[,x] == 1,na.rm=T)/sum(!is.na(genon[,x])))
    toInfer <- (MAF > filter$MAF) & (miss < filter$MISS) & is.na(config)
    
    seg_Infer <- sapply(1:nSnps, function(x){
      if(!toInfer[x])
        return(NA)
      else{
        d = depth_Ref[,x] + depth_Alt[,x]
        g = genon[,x]
        K = sum(1/2^(d[which(d != 0)])*0.5)/sum(d != 0)
        nAA = sum(g==2, na.rm=T)
        nAB = sum(g==1, na.rm=T)
        nBB = sum(g==0, na.rm=T)
        ## check that there are sufficient data to perform the chisq test
        if(sum(nAA+nAB+nBB)/length(g) <= (1-filter$MISS))
          return(NA)
        ## compute chiseq test for both loci types
        exp_prob_BI <- c(0.25 + K,0.5 - 2*K, 0.25 + K)
        exp_prob_SI <- c(K, 0.5 - 2*K, 0.5 + K)
        ctest_BI <- chisq.test(c(nBB,nAB,nAA), p = exp_prob_BI)
        ctest_SI_1 <- chisq.test(c(nBB,nAB,nAA), p = exp_prob_SI)
        ctest_SI_2 <- chisq.test(c(nBB,nAB,nAA), p = rev(exp_prob_SI))
        ## do tests to see if we can infer type
        if( ctest_BI$p.value > filter$PVALUE & ctest_SI_1$p.value < filter$PVALUE & ctest_SI_2$p.value < filter$PVALUE )
          return(1)
        else if ( ctest_BI$p.value < filter$PVALUE & ctest_SI_1$p.value > filter$PVALUE & ctest_SI_2$p.value < filter$PVALUE )
          return(4)
        else if ( ctest_BI$p.value < filter$PVALUE & ctest_SI_1$p.value < filter$PVALUE & ctest_SI_2$p.value > filter$PVALUE )
          return(5)
        else
          return(NA)
      }
    },simplify = T)
    
    indx[[fam]] <- (MAF > filter$MAF) & (miss < filter$MISS) & ( !is.na(config) | !is.na(seg_Infer) )
    
    config[!indx[[fam]]] <- seg_Infer[!indx[[fam]]] <- NA
    
    ## Determine the segregation groups
    config_all[[fam]] <- config
    config_infer_all[[fam]] <- seg_Infer
    
    nSnps_all[[fam]] <- sum(indx[[fam]])
    nInd_all[[fam]] <- nInd
  }
  
  ## Find all the SNPs to keep and subset the global variables
  indx_all <- do.call("rbind",indx)
  indx_all <- apply(indx_all, 2, any)
  
  genon_all     <- genon[,indx_all]
  depth_Ref_all <- depth_Ref[,indx_all]
  depth_Alt_all <- depth_Alt[,indx_all]
  chrom_all     <- R6obj$.__enclos_env__$private$chrom[indx_all]
  pos_all       <- R6obj$.__enclos_env__$private$pos[indx_all]
  SNP_Names     <- R6obj$.__enclos_env__$private$SNP_Names[indx_all]
  
  group <- group_infer <- vector(mode="list", length=noFam)
  for(fam in 1:noFam){
    group[[fam]]$BI <- which(config_all[[fam]][indx_all] == 1)
    group[[fam]]$PI <- which(config_all[[fam]][indx_all] %in% c(2,3))
    group[[fam]]$MI <- which(config_all[[fam]][indx_all] %in% c(4,5))
    
    group_infer[[fam]]$BI <- which(config_infer_all[[fam]][indx_all] == 1) 
    group_infer[[fam]]$SI <- which(config_infer_all[[fam]][indx_all] %in% c(4,5))
    
    config_all[[fam]] <- config_all[[fam]][indx_all]
    config_infer_all[[fam]] <- config_infer_all[[fam]][indx_all]
    
    cat("-------------\n")
    cat("Family ",names(famInfo)[fam]," Summary:\n\n",sep="")
    cat("Number of SNPs remaining after filtering:",nSnps_all[[fam]],"\n")
    cat("Number of progeny:", nInd_all[[fam]],"\n")
    cat("Number of SNPs with correct segregation type:", sum(!is.na(config_all[[fam]])) ,"\n")
    cat("Both-informative (BI):", length(group[[fam]]$BI),"\n")
    cat("Maternal-informative (MI):", length(group[[fam]]$MI),"\n")
    cat("Paternal-informative (PI):", length(group[[fam]]$PI),"\n")
    cat("Number of SNPs with inferred segregation type:", sum(!is.na(config_infer_all[[fam]])),"\n")
    cat("Both-informative (BI):", length(group_infer[[fam]]$BI),"\n")
    cat("Maternal/Paternal-informative (MI or PI):", length(group_infer[[fam]]$SI),"\n")
  }
  
  ## Update the R6 object and return it
  R6obj$.__enclos_env__$private$updatePrivate(list(
    genon = genon, depth_Ref = depth_Ref_all, depth_Alt = depth_Alt_all, chrom = chrom_all, pos = pos_all,
    group = group, group_infer = group_infer, config = config_all, config_infer = config_infer_all,
    nInd = nInd_all, nSnps = nSnps_all, noFam = noFam, indID = indID_all, SNP_Names = SNP_Names)
  )
  
  return(R6obj)
}
