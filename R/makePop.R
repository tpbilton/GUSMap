

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
makePop.FS <- function(R6obj, famInfo, filter=list(MAF=0.05, MISS=0.2, BIN=0, DEPTH=6, PVALUE=0.05), inferSNPs = FALSE, perInfFam=1){
  
  ## Do some checks
  if(is.null(filter$MAF) || filter$MAF<0 || filter$MAF>1 || !is.numeric(filter$MAF))
    stop("Minor allele frequency filter has not be specifies or is invalid.")
  if(perInfFam <=0.5 || perInfFam > 1)
    stop("The of the percentage of families which are informative for each SNP must greater than 50% and leass than equal to 100%.")
  
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
  
  genon_all <- depth_Ref_all <- depth_Alt_all <- vector(mode="list", length=noFam)
  
  ## extract the data and format correct for each family.
  for(fam in 1:noFam){
    cat("Processing Family ",names(famInfo)[fam],".\n\n",sep="")
    mum <- famInfo[[fam]]$parents$Mother
    dad <- famInfo[[fam]]$parents$Father
    patgrandmum <- famInfo[[fam]]$grandparents$paternalGrandMother
    patgranddad <- famInfo[[fam]]$grandparents$paternalGrandFather
    matgrandmum <- famInfo[[fam]]$grandparents$maternalGrandMother
    matgranddad <- famInfo[[fam]]$grandparents$maternalGrandFather
    ## index the parents
    mumIndx <- which(indID %in% mum)
    if(length(mumIndx) == 0)
      stop(paste0("Mother ID not found family ",fam,"."))
    dadIndx <- which(indID %in% dad)
    if(length(dadIndx) == 0)
      stop(paste0("Father ID not found family ",fam,"."))
    ## index the grandparents
    patgrandparents <- matgrandparents <- FALSE
    patgrandmumIndx <- which(indID %in% patgrandmum)
    patgranddadIndx <- which(indID %in% patgranddad)
    matgrandmumIndx <- which(indID %in% matgrandmum)
    matgranddadIndx <- which(indID %in% matgranddad)
    if(!is.null(patgrandmumIndx) && !is.null(patgranddadIndx))
      patgrandparents <- TRUE
    if(!is.null(matgrandmumIndx) && !is.null(matgranddadIndx))
      matgrandparents <- TRUE
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
    
    if(patgrandparents){
      genon_patgrandmum <- matrix(R6obj$.__enclos_env__$private$genon[patgrandmumIndx,], nrow=length(patgrandmumIndx), ncol=nSnps) 
      depth_patgrandmum <- matrix(R6obj$.__enclos_env__$private$depth_Ref[patgrandmumIndx,] +
                            R6obj$.__enclos_env__$private$depth_Alt[patgrandmumIndx,], nrow=length(patgrandmumIndx), ncol=nSnps)
      genon_patgranddad <- matrix(R6obj$.__enclos_env__$private$genon[patgranddadIndx,], nrow=length(patgranddadIndx), ncol=nSnps) 
      depth_patgranddad <- matrix(R6obj$.__enclos_env__$private$depth_Ref[patgranddadIndx,] +
                                    R6obj$.__enclos_env__$private$depth_Alt[patgranddadIndx,], nrow=length(patgranddadIndx), ncol=nSnps)
    }
    if(matgrandparents){
      genon_matgrandmum <- matrix(R6obj$.__enclos_env__$private$genon[matgrandmumIndx,], nrow=length(matgrandmumIndx), ncol=nSnps) 
      depth_matgrandmum <- matrix(R6obj$.__enclos_env__$private$depth_Ref[matgrandmumIndx,] +
                                    R6obj$.__enclos_env__$private$depth_Alt[matgrandmumIndx,], nrow=length(matgrandmumIndx), ncol=nSnps)
      genon_matgranddad <- matrix(R6obj$.__enclos_env__$private$genon[matgranddadIndx,], nrow=length(matgranddadIndx), ncol=nSnps) 
      depth_matgranddad <- matrix(R6obj$.__enclos_env__$private$depth_Ref[matgranddadIndx,] +
                                    R6obj$.__enclos_env__$private$depth_Alt[matgranddadIndx,], nrow=length(matgranddadIndx), ncol=nSnps)
    }
    
    parHap_pat <- sapply(1:nSnps,function(x){
      x_p = genon_dad[,x]; d_p = depth_dad[,x]
      if(any(x_p==1,na.rm=T))
        return("AB")
      else if(sum(d_p) > filter$DEPTH){
        if(x_p==2)
          return("AA")
        else if(x_p==0)
          return("BB")
        else if(patgrandparents){
          if(sum(depth_patgranddad[,x])>filter$DEPTH && sum(depth_patgrandmum)>filter$DEPTH){
            x_gp = genon_patgranddad[,x]; x_gm = genon_patgrandmum[,x]
            if((x_gp == 2 & x_gm == 0) || (x_gp == 0 & x_gm == 2))
              return("AB")
            else if(x_gp == 2 & x_gm == 2)
              return("AA")
            else if(x_gp == 0 & x_gm == 0)
              return("BB")
          }
        }
        else
          return(NA)
      }
      else
        return(NA)
    })
    
    parHap_mat <- sapply(1:nSnps,function(x){
      x_m = genon_mum[,x]; d_m = depth_mum[,x]
      if(any(x_m==1,na.rm=T))
        return("AB")
      else if(sum(d_m) > filter$DEPTH){
        if(x_m==2)
          return("AA")
        else if(x_m==0)
          return("BB")
        else if(matgrandparents){
          if(sum(depth_matgranddad[,x])>filter$DEPTH && sum(depth_matgrandmum)>filter$DEPTH){
            x_gp = genon_matgranddad[,x]; x_gm = genon_matgrandmum[,x]
            if((x_gp == 2 & x_gm == 0) || (x_gp == 0 & x_gm == 2))
              return("AB")
            else if(x_gp == 2 & x_gm == 2)
              return("AA")
            else if(x_gp == 0 & x_gm == 0)
              return("BB")
          }
        }
        else
          return(NA)
      }
      else
        return(NA)
    })
    
    config <- rep(NA,nSnps)
    config[which(parHap_pat == "AB" & parHap_mat == "AB")] <- 1
    config[which(parHap_pat == "AB" & parHap_mat == "AA")] <- 2
    config[which(parHap_pat == "AB" & parHap_mat == "BB")] <- 3
    config[which(parHap_pat == "AA" & parHap_mat == "AB")] <- 4
    config[which(parHap_pat == "AA" & parHap_mat == "BB")] <- 5
    
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
          ctest <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = exp_prob))
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
        else if(config[x] %in% c(2,4)){
          exp_prob <- c(K, 0.5 - 2*K, 0.5 + K)
          ctest <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = exp_prob))
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
        else if(config[x] %in% c(3,5)){
          exp_prob <- c(0.5 + K, 0.5 - 2*K, K)
          ctest <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = exp_prob))
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
    if(inferSNPs){
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
          ctest_BI <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = exp_prob_BI))
          ctest_SI_1 <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = exp_prob_SI))
          ctest_SI_2 <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = rev(exp_prob_SI)))
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
    }
    
    chrom <- R6obj$.__enclos_env__$private$chrom
    pos <- R6obj$.__enclos_env__$private$pos
    ## Extract one SNP from each read.
    if(filter$BIN > 0){
      oneSNP <- rep(FALSE,nSnps)
      oneSNP[unlist(sapply(unique(chrom), function(x){
        ind <- which(chrom == x)
        g1_diff <- diff(pos[ind])
        SNP_bin <- c(0,cumsum(g1_diff > filter$BIN)) + 1
        set.seed(58473+as.numeric(which(x==chrom))[1])
        keepPos <- sapply(unique(SNP_bin), function(y) {
          ind2 <- which(SNP_bin == y)
          if(length(ind2) > 1)
            return(sample(ind2,size=1))
          else if(length(ind2) == 1)
            return(ind2)
        })
        return(ind[keepPos])
      },USE.NAMES = F ))] <- TRUE
    }
    else 
      oneSNP <- rep(TRUE,nSnps)
    
    if(inferSNPs){
      indx[[fam]] <- (MAF > filter$MAF) & (miss < filter$MISS) & ( !is.na(config) | !is.na(seg_Infer) ) & oneSNP
      config[!indx[[fam]]] <- seg_Infer[!indx[[fam]]] <- NA
    }
    else
      indx[[fam]] <- (MAF > filter$MAF) & (miss < filter$MISS) & (!is.na(config)) & oneSNP
  
    
    ## Determine the segregation groups
    config_all[[fam]] <- config
    if(inferSNPs)
      config_infer_all[[fam]] <- seg_Infer

    nSnps_all[[fam]] <- sum(indx[[fam]])
    nInd_all[[fam]] <- nInd
    
    genon_all[[fam]] <- genon
    depth_Ref_all[[fam]] <- depth_Ref
    depth_Alt_all[[fam]] <- depth_Alt
  }
  
  if(noFam == 1){
    fam = 1
    
    ## Find all the SNPs to keep and subset the global variables
    indx_all <- indx[[fam]]
    #indx_all <- do.call("rbind",indx)
    #indx_all <- apply(indx_all, 2, any)
    
    genon_all[[fam]]     <- genon[,indx_all]
    depth_Ref_all[[fam]] <- depth_Ref[,indx_all]
    depth_Alt_all[[fam]] <- depth_Alt[,indx_all]
    chrom_all            <- R6obj$.__enclos_env__$private$chrom[indx_all]
    pos_all              <- R6obj$.__enclos_env__$private$pos[indx_all]
    SNP_Names            <- R6obj$.__enclos_env__$private$SNP_Names[indx_all]
    
    #group <- group_infer <- vector(mode="list", length=noFam)
    group <- group_infer <- list()
    group$BI <- which(config_all[[fam]][indx_all] == 1)
    group$PI <- which(config_all[[fam]][indx_all] %in% c(2,3))
    group$MI <- which(config_all[[fam]][indx_all] %in% c(4,5))
      
    group_infer$BI <- which(config_infer_all[[fam]][indx_all] == 1) 
    group_infer$SI <- which(config_infer_all[[fam]][indx_all] %in% c(4,5))
      
    config_all[[fam]] <- config_all[[fam]][indx_all]
    if(inferSNPs)
      config_infer_all[[fam]] <- config_infer_all[[fam]][indx_all]
      
    cat("-------------\n")
    cat("Summary:\n\n")
    cat("Number of SNPs remaining after filtering:",sum(indx_all),"\n")
    if(inferSNPs)
      cat("Number of SNPs with correct segregation type:", sum(!is.na(config_all[[fam]])) ,"\n")
    cat("Both-informative (BI):", length(group$BI),"\n")
    cat("Maternal-informative (MI):", length(group$MI),"\n")
    cat("Paternal-informative (PI):", length(group$PI),"\n")
    if(inferSNPs){
      cat("Number of SNPs with inferred segregation type:", sum(!is.na(config_infer_all[[fam]])),"\n")
      cat("Both-informative (BI):", length(group_infer$BI),"\n")
      cat("Maternal/Paternal-informative (MI or PI):", length(group_infer$SI),"\n")
    }
    cat("Number of progeny:", nInd_all[[fam]],"\n")
  }
  else{
    noInfoFam <- ceiling(perInfFam*noFam)
    group <- group.temp <- list()
    ## Work out the SNP groupings
    # BI
    tabInf_BI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(1) & indx[[y]]))))
    group.temp$BI <- as.numeric(names(tabInf_BI)[which(tabInf_BI >= noInfoFam)])
    # MI 
    tabInf_MI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(1,4,5) & indx[[y]]))))
    group.temp$MI <- setdiff(as.numeric(names(tabInf_MI)[which(tabInf_MI >= noInfoFam)]), group.temp$BI)
    # PI 
    tabInf_PI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(1,2,3) & indx[[y]]))))
    group.temp$PI <- setdiff(as.numeric(names(tabInf_PI)[which(tabInf_PI >= noInfoFam)]), group.temp$BI)
    
    ## Work out which SNPs to keep
    indx_all <- logical(nSnps)
    indx_all[unlist(group.temp)] <- TRUE
    
    add.bi <- group.temp$PI[which(group.temp$PI %in% group.temp$MI)]
    group.temp$BI <- sort(c(group.temp$BI,add.bi))
    group.temp$MI <- setdiff(group.temp$MI,group.temp$BI)
    group.temp$PI <- setdiff(group.temp$PI,group.temp$BI)
    
    group$BI <- which(sort(unique(unlist(group.temp))) %in% group.temp$BI)
    group$MI <- which(sort(unique(unlist(group.temp))) %in% group.temp$MI)
    group$PI <- which(sort(unique(unlist(group.temp))) %in% group.temp$PI)
    
    ## Subset the data
    for(fam in 1:noFam){
      genon_all[[fam]]     <- genon_all[[fam]][,indx_all]
      depth_Ref_all[[fam]] <- depth_Ref_all[[fam]][,indx_all]
      depth_Alt_all[[fam]] <- depth_Alt_all[[fam]][,indx_all]
      config_all[[fam]] <- config_all[[fam]][indx_all]
    }
    chrom_all     <- R6obj$.__enclos_env__$private$chrom[indx_all]
    pos_all       <- R6obj$.__enclos_env__$private$pos[indx_all]
    SNP_Names     <- R6obj$.__enclos_env__$private$SNP_Names[indx_all]
    
    cat("-------------\n")
    cat("Summary:\n\n",sep="")
    cat("Number of SNPs remaining after filtering:",sum(indx_all),"\n")
    cat("Both-informative (BI):", length(group$BI),"\n")
    cat("Maternal-informative (MI):", length(group$MI),"\n")
    cat("Paternal-informative (PI):", length(group$PI),"\n")
    cat("Number of progeny in Family...\n")
    for(fam in 1:noFam)
      cat(names(famInfo)[fam],":", nInd_all[[fam]],"\n")
    cat("\n")
  
    group_infer <- NULL
    config_infer_all <- NULL
  }
  
  ## Update the R6 object and return it
  R6obj$.__enclos_env__$private$updatePrivate(list(
    genon = genon_all, depth_Ref = depth_Ref_all, depth_Alt = depth_Alt_all, chrom = chrom_all, pos = pos_all,
    group = group, group_infer = group_infer, config = config_all, config_infer = config_infer_all,
    nInd = nInd_all, nSnps = nSnps_all, noFam = noFam, indID = indID_all, SNP_Names = SNP_Names)
  )
  
  return(R6obj)
}
