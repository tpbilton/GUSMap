#' @export FS

### R6 class for creating a data format for full-sib families
FS <- R6Class("FS",
              inherit = RA,
              public = list(
                ## initialize function
                initialize = function(R6obj){
                  private$masked    <- rep(FALSE, R6obj$.__enclos_env__$private$nSnps)
                  private$famInfo   <- R6obj$.__enclos_env__$private$famInfo
                  private$genon     <- R6obj$.__enclos_env__$private$genon
                  private$ref       <- R6obj$.__enclos_env__$private$ref
                  private$alt       <- R6obj$.__enclos_env__$private$alt
                  private$chrom     <- R6obj$.__enclos_env__$private$chrom
                  private$pos       <- R6obj$.__enclos_env__$private$pos
                  private$SNP_Names <- R6obj$.__enclos_env__$private$SNP_Names
                  private$indID     <- R6obj$.__enclos_env__$private$indID
                  private$nSnps     <- R6obj$.__enclos_env__$private$nSnps
                  private$nInd      <- R6obj$.__enclos_env__$private$nInd
                  private$gform     <- R6obj$.__enclos_env__$private$gform
                  private$AFrq      <- R6obj$.__enclos_env__$private$AFrq
                  private$infilename<- R6obj$.__enclos_env__$private$infilename
                  private$para      <- list(OPGP=list(),rf_p=list(),rf_m=list(),ep=list(),loglik=list())
                },
                ### print methods
                 print = function(...){
                   if(noFam == 1){
                     cat("Single Family Linkage analysis:\n\n")
                     cat("Data Summary:\n")
                     cat("Data file:\t",private$infilename,"\n")
                     temp <- private$ref[[1]] + private$alt[[1]]
                     cat("Mean Depth:\t", mean(temp),"\n")
                     cat("Mean Call Rate:\t",sum(temp!=0)/length(temp),"\n")
                     cat("Number of ...\n")
                     cat("  Progeny:\t",unlist(private$nInd),"\n")
                     cat("  MI SNPs:\t",length(private$group$MI),"\n")
                     cat("  PI SNPs:\t",length(private$group$PI),"\n")
                     cat("  BI SNPs:\t",length(private$group$BI),"\n")
                     cat("  Total SNPs:\t",length(unlist(private$group)),"\n\n")
                     ## linkage group information (if any)
                     if(is.null(private$LG)){
                       cat("Linkage Group Summary:")
                       MI <- unlist(lapply(private$LG_mat, function(x) sum(x %in% private$group$MI)))
                       PI <- unlist(lapply(private$LG_pat, function(x) sum(x %in% private$group$PI)))
                       tab <- cbind(LG=1:(length(MI)+length(PI)),MI=c(MI,rep(0,length(PI))),PI=c(rep(0,length(MI)),PI))
                       prmatrix(tab, rowlab = rep("",nrow(tab)))
                     }
                     else{
                       cat("Linkage Group Summary:")
                       MI <- unlist(lapply(private$LG, function(x) sum(x %in% private$group$MI)))
                       PI <- unlist(lapply(private$LG, function(x) sum(x %in% private$group$PI)))
                       BI <- unlist(lapply(private$LG, function(x) sum(x %in% private$group$BI)))
                       TOTAL <- unlist(lapply(private$LG, length))
                       tab <- cbind(LG=1:(length(private$LG)),MI,PI,BI,TOTAL)
                       if(!is.null(private$para)){
                        DIST_MAT <- extendVec(unlist(lapply(private$para$rf_m, function(x) round(sum(mfun(x, fun="haldane", centiM = T)),2)), nrow(tab)))
                        DIST_PAT <- extendVec(unlist(lapply(private$para$rf_p, function(x) round(sum(mfun(x, fun="haldane", centiM = T)),2)), nrow(tab)))
                        ERR <- extendVec(round(unlist(private$para$ep),5), nrow(tab))
                        tab <- cbind(tab,DIST_MAT,DIST_PAT,ERR)
                        prmatrix(tab, rowlab = rep("",nrow(tab)))
                       }
                     }
                   }
                   else
                     cat("Not yet implemented")
                   invisible(self)
                 },
                #############################################################
                ## Function for removing SNPs from the linkage groups
                removeSNP = function(indx){
                  ## some checks
                  if( !is.vector(indx) || !is.numeric(indx) || !all(indx == round(indx)) || any(indx < 1) || any(indx > private$nSnps) )
                    stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                  indx <- unique(indx) ## make sure no double ups in SNPs
                  ## remove SNP from the maternal linkage groups
                  for(lg in 1:length(private$LG_mat)){
                    if(any(private$LG_mat[[lg]] %in% indx))
                      private$LG_mat[[lg]] <- private$LG_mat[[lg]][-which(private$LG_mat[[lg]] %in% indx)]
                  }
                  ## remove SNP from the paternal linkage groups
                  for(lg in 1:length(private$LG_pat)){
                    if(any(private$LG_pat[[lg]] %in% indx))
                      private$LG_pat[[lg]] <- private$LG_pat[[lg]][-which(private$LG_pat[[lg]] %in% indx)]
                  }
                  return(invisible())
                },
                ## function for masking SNPs
                maskSNP = function(indx){
                  ## Inout checks
                  if( !is.vector(indx) || !is.numeric(indx) || !all(indx == round(indx)) || any(indx < 1) || any(indx > private$nSnps) )
                    stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                  ## mask SNPs
                  private$masked[indx] <- TRUE
                  return(self)
                },
                ## function for unmasking SNPs
                unMaskSNP = function(indx){
                  if( !is.vector(indx) || !is.numeric(indx) || !all(indx == round(indx)) || any(indx < 1) || any(indx > private$nSnps) )
                    stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                  private$masked[indx] <- FALSE
                  return(self)
                },
                #####################################################################
                ## Function for computing the 2-point rf estimates
                rf_2pt = function(nClust=4){
                  ## do some checks
                  if(!is.numeric(nClust) || !is.vector(nClust) || length(nClust)!=1 || nClust < 0 || is.infinite(nClust) )
                    stop("Number of clusters for parallelization needs to be a positive finite integer number")
                  
                  ## If there is only one family
                  if(private$noFam == 1)
                    mat <- rf_2pt_single(private$ref[[1]], private$alt[[1]],
                                         private$config[[1]], private$config_infer[[1]],
                                         private$group, private$group_infer,
                                         nClust, private$nInd)
                  else
                    mat <- rf_2pt_multi(private$ref, private$alt,
                                        private$config,private$group, nClust, private$noFam)
                  ## Save the results to the object
                  private$rf <- mat$rf
                  private$LOD <- mat$LOD
                  return(invisible())
                },
                ## Function for creating linkage groups
                createLG = function(parent, LODthres=10, nComp=10){
                  ## Do some checks
                  if(is.null(private$rf) || is.null(private$LOD))
                    stop("Recombination fractions and LOD scores have not been computed.\nUse rf_2pt() to compute the recombination fractions and LOD scores.")
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop("parent argument is not a string of length one or is invalid:
Please select one of the following:
   maternal: Only MI SNPs
   paternal: Only PI SNPs
   both:     MI and PI SNPs")
                  if(!is.numeric(LODthres) || !is.vector(LODthres) || length(LODthres) != 1 || LODthres < 0 || !is.finite(LODthres))
                    stop("The LOD threshold (argument 2) needs to be a finite numeric number.")
                  if(!is.numeric(nComp) || !is.vector(nComp) || length(nComp) != 1 || nComp < 0 || !is.finite(nComp) ||
                     round(nComp) != nComp)
                    stop("The number of comparsion points (argument 3) needs to be a finite integer number.")
                  
                  ## Create the groups
                  if(parent == "maternal" || parent == "both")
                    private$LG_mat <- createLG(private$group, private$LOD, "maternal", LODthres, nComp, masked=private$masked)
                  if(parent == "paternal" || parent == "both")
                    private$LG_pat <- createLG(private$group, private$LOD, "paternal", LODthres, nComp, masked=private$masked)
                  return(invisible())
                },
                ## function for adding the unmapped (or inferred SNPs) to the linkage groups
                # addSNPs = function(LODthres=10, nComp=10){
                # 
                #   if(length(self$LG_mat) == 0 || length(self$LG_pat) == 0)
                #     stop("There are no linkage groups. Please use 'createLG' function to create the linkage groups first")
                #   if(!is.numeric(LODthres) || !is.vector(LODthres) || length(LODthres) != 1 || LODthres < 0 || !is.finite(LODthres))
                #     stop("The LOD threshold (argument 1) needs to be a finite positive numeric number.")
                #   if(!is.numeric(nComp) || !is.vector(nComp) || length(nComp) != 1 || nComp < 0 || !is.finite(nComp) ||
                #      round(nComp) != nComp)
                #     stop("The number of comparsion points (argument 2) needs to be a finite positive integer number.")
                # 
                #   ## Find the unmapped loci
                #   
                #   unmapped <- sort(unlist(c(which(!(private$group$MI %in% unlist(c(self$LG_mat,self$LG_pat)))),
                #                             which(!(private$group$PI %in% unlist(c(self$LG_mat,self$LG_pat)))),
                #                             private$group_infer$SI)))
                #   if(length(unmapped) == 0)
                #     stop("There are no SNPs remaining that are unmapped")
                # 
                #   ## new LGs for adding BI snps to
                #   newLGlist <- c(self$LG_mat, self$LG_pat)
                #   ## count the number of LGs
                #   nLG <- length(newLGlist)
                #   
                #   ## Run algorithm for generating the linkage groups
                #   noneMapped = FALSE
                #   count = 0
                #   while(!noneMapped){
                #     noneMapped = TRUE
                #     ## check that there are still SNPs remaining that need to be mapped
                #     if(length(unmapped) == 0)
                #       next
                #     ## run the algorithm to map the SNPs
                #     else{
                #       for(snp in unmapped){
                #         LODvalue = numeric(nLG)
                #         for(lg in 1:nLG)
                #           LODvalue[lg] <- mean(sort(private$LOD[snp,newLGlist[[lg]]],decreasing=T)[1:nComp],na.rm=T)
                #         if(max(LODvalue) >= LODthres & sort(LODvalue, decreasing = T)[2] < LODthres){
                #           count = count + 1
                #           newLG <- which.max(LODvalue)
                #           newLGlist[[newLG]] <- c(newLGlist[[newLG]], snp)
                #           unmapped <- unmapped[-which(unmapped == snp)]
                #           noneMapped = FALSE
                #         }
                #       }
                #     }
                #   }
                #   ## Check to see if any SNPs were added
                #   if(count == 0){
                #     stop("No SNPs were added to the Linkage Groups")
                #   } else{ ## else set the new LGs
                #     self$LG_mat_temp <- newLGlist[1:length(self$LG_mat)]
                #     self$LG_pat_temp <- newLGlist[length(self$LG_mat) + 1:length(self$LG_pat)]
                #   }
                #   return(invisible())
                # },
                ## Function for adding the informative SNPs to the LGs
                addBIsnps = function(parent = "maternal", LODthres=10, nComp=10){

                  ## Do some checks
                  if(is.null(private$rf) || is.null(private$LOD))
                    stop("Recombination fractions and LOD scores have not been computed.\nUse rf_2pt() to compute the recombination fractions and LOD scores.")
                  if((parent == "maternal" || parent == "both") && is.null(private$LG_mat))
                    stop("There are no maternal linkage groups. Use createLG(parent='maternal') to form maternal linkage groups.")
                  if((parent == "paternal" || parent == "both") && is.null(private$LG_pat))
                    stop("There are no maternal linkage groups. Use createLG(parent='maternal') to form maternal linkage groups.")
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop("parent argument is not a string of length one or is invalid:
Please select one of the following:
   maternal: Add BI SNPs to MI LGs
   paternal: Add BI SNPs to PI LGs
   both:     Add BI SNPs to both MI and PI LGs")
                  if(parent == "maternal" && length(private$LG_mat) == 0)
                    stop("There are no maternal linkage groups. Please use 'createLG' function to create the linkage groups first")
                  else if(parent == "paternal" && length(private$LG_pat) == 0)
                    stop("There are no paternal linkage groups. Please use 'createLG' function to create the linkage groups first")
                  if(!is.numeric(LODthres) || !is.vector(LODthres) || length(LODthres) != 1 || LODthres < 0 || !is.finite(LODthres))
                    stop("The LOD threshold (argument 1) needs to be a finite positive numeric number.")
                  if(!is.numeric(nComp) || !is.vector(nComp) || length(nComp) != 1 || nComp < 0 || !is.finite(nComp) ||
                     round(nComp) != nComp)
                    stop("The number of comparsion points (argument 2) needs to be a finite positive integer number.")
                  
                  ## new LGs for adding BI snps to
                  #if(parent == "maternal")
                  #  newLGlist <- private$LG_mat
                  #else
                  #  newLGlist <- private$LG_pat
                  ## count the number of LGs
                  #nLG <- length(newLGlist)
                  
                  ## Find the unmapped loci
                  unmapped <- sort(unlist(private$group$BI,private$group_infer$BI))
                  ## Remove masked SNPs
                  unmapped <- unmapped[which(!private$masked[unmapped])]
                  ## check that there are SNPs to map
                  if(length(unmapped) == 0)
                    stop("There are no SNPs remaining that are unmapped")
                  
                  if(!is.null(private$LG_mat))
                    newLGlist_mat <- private$mapBISnps(unmapped, parent="maternal", LODthres=LODthres, nComp=nComp)
                  
                  if(!is.null(private$LG_pat))
                    newLGlist_pat <- private$mapBISnps(unmapped, parent="paternal", LODthres=LODthres, nComp=nComp)
                  
                  ## work out which LGs go together
                  indx_mat <- sapply(unmapped, function(x){ which(unlist(lapply(newLGlist_mat, function(y) any(x == y))))})
                  indx_pat <- sapply(unmapped, function(x){ which(unlist(lapply(newLGlist_pat, function(y) any(x == y))))})
                  indx_new <- clue::solve_LSAP(table(indx_mat,indx_pat), maximum=T)
                  indx_mat <- car::recode(indx_mat, paste(paste(1:5,indx_new,sep="="),collapse = ";"))
                  
                  ## remove any BI SNPs that are problematic
                  
                  ## Run algorithm for generating the linkage groups
                  # noneMapped = FALSE
                  # count = 0
                  # while(!noneMapped){
                  #   noneMapped = TRUE
                  #   ## check that there are still SNPs remaining that need to be mapped
                  #   if(length(unmapped) == 0)
                  #     next
                  #   ## run the algorithm to map the SNPs
                  #   else{
                  #     for(snp in unmapped){
                  #       LODvalue = numeric(nLG)
                  #       for(lg in 1:nLG)
                  #         LODvalue[lg] <- mean(sort(private$LOD[snp,newLGlist[[lg]]],decreasing=T)[1:nComp],na.rm=T)
                  #       if(max(LODvalue) >= LODthres & sort(LODvalue, decreasing = T)[2] < LODthres){
                  #         count = count + 1
                  #         newLG <- which.max(LODvalue)
                  #         newLGlist[[newLG]] <- c(newLGlist[[newLG]], snp)
                  #         unmapped <- unmapped[-which(unmapped == snp)]
                  #         noneMapped = FALSE
                  #       }
                  #     }
                  #   }
                  #}
                  if(count == 0){
                    stop("No SNPs were added to the Linkage Groups")
                  } else{
                    if(parent == "maternal")
                      private$LG_mat_temp <- newLGlist_mat
                    else
                      private$LG_pat_temp <- newLGlist_pat
                  }
                  return(invisible())
                },
                ## Function for ordering linkage groups
                orderLG = function(chrom = NULL, mapfun = "morgan", weight="LOD2", ndim=2, spar=NULL){
                  ## do some checks
                  if((length(self$LG_mat) != length(self$LG_pat)))
                    stop("The number of linkage groups for the two parental meiosis are not the same")
                  if(is.null(chrom))
                    chrom <- 1:length(self$LG_mat)
                  else if(!is.vector(chrom) && !is.numeric(chrom) && any(chrom < 0) && chrom > length(self$LG_mat))
                    stop("invalid chromosome input")
                  if(!(mapfun %in% c("morgan","haldane","kosambi")))
                    stop("Unknown mapping function")
                  if(!(weight %in% c("LOD","LOD2","none")))
                     stop("Unknown weighting function")
                  nChr = length(chrom)
                  LG <- vector(mode='list', length=nChr)
                  ## order the linkage groups
                  for(chr in chrom){
                    ind <- unique(c(self$LG_mat[[chr]], self$LG_pat[[chr]]))
                    ind_mi <- which(ind %in% private$group$MI)
                    ind_pi <- which(ind %in% private$group$PI)
                    ## set up the weighting matrix
                    if(weight == "LOD")
                      wmat <- matrix(private$LOD,nrow=length(ind), ncol=length(ind))
                    else if (weight == "LOD2")
                      wmat <- matrix((private$LOD)^2,nrow=length(ind), ncol=length(ind))
                    else if (weight == "none")
                      wmat <- matrix(1,nrow=length(ind), ncol=length(ind))
                    ## set the weights of the PI and MI combinations to zero.
                    wmat[ind_mi, ind_pi] <- 0
                    wmat[ind_pi, ind_mi] <- 0
                    ## set of the distance matrix
                    dmat <- mfun(private$rf[ind,ind], fun = mapfun)
                    dmat[which(is.infinite(dmat))] <- 18.3684
                    ## unconstrainted MDS
                    MDS <- smacof::smacofSym(delta=dmat, ndim=ndim, weightmat = wmat, itmax = 1e+05)
                    ## principal curve
                    pcurve <- princurve::principal_curve(MDS$conf, maxit = 150, spar = spar)
                    ## plot the results
                    par(mfrow = c(1,2))
                    graphics::plot(MDS$conf[,1],MDS$conf[,2], ylab="Dimension 2", xlab="Dimension 1", type="n")
                    text(MDS$conf, labels=ind)
                    lines(pcurve)
                    graphics::image(private$rf[ind,ind][pcurve$ord,pcurve$ord], axes=F)
                    ## Set the new order
                    self$LG[[chr]] <- ind[pcurve$ord]
                  }
                  return(invisible)
                },
                ## Function for setting temp LGs to new LGs
                setLG = function(parent="both"){
                  if(parent == "both"){
                    private$LG_mat <- private$LG_mat_temp
                    private$LG_mat_temp <- NULL
                    private$LG_pat <- private$LG_pat_temp
                    private$LG_pat_temp <- NULL
                    cat("Temporary linkage groups now set as new linkage groups\n")
                  }
                  else if(parent == "maternal"){
                    private$LG_mat <- private$LG_mat_temp
                    private$LG_mat_temp <- NULL
                    cat("Temporary maternal linkage groups now set as new linkage groups\n")
                  }
                  else if(parent == "paternal"){
                    private$LG_pat <- private$LG_pat_temp
                    private$LG_pat_temp <- NULL
                    cat("Temporary paternal linkage groups now set as new linkage groups\n")
                  }
                  return(invisible())
                },
                ## Function for plotting linkage groups
                plotLG = function(mat=c("rf"), LG, filename=NULL, names=NULL, chrS=2, lmai=2, chrom=T){
                  
                  if(mat == "rf")
                    plotLG(mat=private$rf, LG=LG, filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=chrom)
                  else if(mat == "LOD")
                    plotLG(mat=private$LOD, LG=LG, filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=chrom)
                  else
                    stop("Matrix to be plotted not found.")
                  return(invisible())
                },
                ## Function for plotting chromosome ordering
                plotChr = function(mat=c("rf"), parent = "maternal", ordering = "original", filename=NULL, chrS=2, lmai=2){
                  if(private$noFam == 1){
                    ## workout the indices for the chromosomes
                    names <- unique(private$chrom)
                    if(parent == "maternal"){
                      if(ordering == "original")
                        LG <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,4,5))), simplify=F)
                      else if(ordering == "LG")
                        LG <- private$LG_mat
                      else if(ordering == "tempLG")
                        LG <- private$LG_mat_temp
                    }
                    if(parent == "paternal"){
                      if(ordering == "original")
                        LG <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,2,3))), simplify=F)
                      else if(ordering == "LG")
                        LG <- private$LG_pat
                      else if(ordering == "tempLG")
                        LG <- private$LG_pat_temp
                    }
                    ## plot the chromsomes rf info
                    if(mat == "rf")
                      plotLG(mat=private$rf, LG=LG, filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=T)
                    else if(mat == "LOD")
                      plotLG(mat=private$LOD, LG=LG, filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=T)
                    else
                      stop("Matrix to be plotted not found.")
                    return(invisible())
                    }
                  else{
                    stop("not implemented yet")
                  }
                },
                ## Function for computing the rf's for each chromosome 
                rf_est = function(chr=NULL, init_r=0.01, ep=0.001, method="optim", sexSpec=F, seqErr=T, mapped=T){
                  ## do some checks
                  if( !is.null(init_r) & !is.numeric(init_r) )
                    stop("Starting values for the recombination fraction needs to be a numeric vector or integer or a NULL object")
                  if( (length(ep) != 1 || !is.numeric(ep) || (ep <= 0 | ep >= 1)) )
                    stop("Value for the error parameters needs to be a single numeric value in the interval (0,1) or a NULL object")
                  ## for existing chromosome orders
                  if(!mapped){
                    if(!is.null(chr) && !is.vector(chr))
                      stop("Chromosomes names must be a character vector.")
                    if(is.null(chr)) # compute distances for all chromosomes
                      chr <- unique(private$chrom)
                    else if(!(any(chr %in% private$chrom[!private$masked])))
                      stop("At least one chromosome not found in the data set.")
                    else
                      chr <- as.character(chr) # ensure that the chromsome names are characters
                    cat("Computing recombination fractions:\n")
                    for(i in chr){
                      cat("Chromosome: ",i,"\n")
                      indx_chr <- which((private$chrom == i) & !private$masked)
                      ref_temp <- lapply(private$ref, function(x) x[,indx_chr])
                      alt_temp <- lapply(private$alt, function(x) x[,indx_chr])
                      ## estimate OPGP's
                      if(!any(names(self$para$OPGP) != i) || sapply(1:private$noFam, function(x)
                        !is.null(self$para$OPGP[i][[x]]) && (length(self$para$OPGP[i][[x]]) != ncol(ref_temp[[x]])))){
                        tempOPGP <- list()
                        for(fam in 1:private$noFam){
                          tempOPGP <- c(tempOPGP,list(as.integer(infer_OPGP_FS(ref_temp[[fam]],alt_temp[[fam]],private$config[[fam]][indx_chr], method="EM"))))                        
                          }
                        self$para$OPGP <- tempOPGP
                      }
                      ## estimate the rf's
                      MLE <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=self$para$OPGP,
                                         sexSpec=sexSpec, seqErr=seqErr, method=method)
                      if(sexSpec){
                        self$para$rf_p[i]   <- list(MLE$rf_p)
                        self$para$rf_m[i]   <- list(MLE$rf_m)
                      } else{
                        self$para$rf_p[i]   <- list(MLE$rf)
                        self$para$rf_m[i]   <- list(MLE$rf)
                      }
                      self$para$ep[i]    <- list(list(MLE$ep))
                      self$para$loglik[i]<- list(MLE$loglik)
                    }
                  }
                  else{
                    if((length(self$LG) == 0) && length(is.null(self$LG) == 0))
                      stop("No linkage groups have been formed.")
                    if(!is.null(chr) && (!is.vector(chr) || chr < 0 || chr > length(self$LG) || round(chr) != chr))
                      stop("Chromosomes input must be an integer vector between 1 and the number of linkage groups")
                    cat("Computing recombination fractions:\n")
                    for(i in chr){
                      cat("Chromosome: ",i,"\n")
                      indx_chr <- self$LG[[chr]]
                      ref_temp <- lapply(private$ref, function(x) x[,indx_chr])
                      alt_temp <- lapply(private$alt, function(x) x[,indx_chr])
                      ## estimate OPGP's
                      if(!any(names(self$para$OPGP) != i) || sapply(1:private$noFam, function(x)
                        !is.null(self$para$OPGP[i][[x]]) && (length(self$para$OPGP[i][[x]]) != ncol(ref_temp[[x]])))){
                        tempOPGP <- list()
                        for(fam in 1:private$noFam){
                          tempOPGP <- c(tempOPGP,list(as.integer(infer_OPGP_FS(ref_temp[[fam]],alt_temp[[fam]],private$config[[fam]][indx_chr], method="EM"))))                        
                        }
                        self$para$OPGP[i] <- tempOPGP
                      }
                      ## estimate the rf's
                      MLE <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=self$para$OPGP[i],
                                       sexSpec=sexSpec, seqErr=seqErr, method=method)
                      if(sexSpec){
                        self$para$rf_p[i]   <- list(MLE$rf_p)
                        self$para$rf_m[i]   <- list(MLE$rf_m)
                      } else{
                        self$para$rf_p[i]   <- list(MLE$rf)
                        self$para$rf_m[i]   <- list(MLE$rf)
                      }
                      self$para$ep[i]     <- list(MLE$ep)
                      self$para$loglik[i] <- list(MLE$loglik)
                    }
                  }
                  return(invisible())
                },
                ## Ratio of alleles for heterozygous genotype calls (observed vs expected)
                # Slightly different than the one in RA class
                cometPlot = function(model="random", alpha=NULL, filename="HeteroPlot", cex=1, maxdepth=500){
                  ref <- alt <- NULL
                  for(fam in 1:private$noFam){
                    ref <- rbind(ref,private$ref[[fam]])
                    alt <- rbind(alt,private$alt[[fam]])
                  }
                  cometPlot(ref, alt, model=model, alpha=alpha, filename=filename, cex=cex, maxdepth=maxdepth)
                },
                ###### Function for extracting variables from the FS object
                extractVar = function(nameList){
                  if(!is.character(nameList) || !is.vector(nameList) || any(is.na(nameList)))
                  res <- NULL
                  for(name in nameList){
                    res <- c(res,list(private[[name]]))
                  }
                  names(res) <- nameList
                  return(res)
                }
                ##############################################################
              ),
              private = list(
                config       = NULL,
                config_infer = NULL,
                group        = NULL,
                group_infer  = NULL,
                masked       = NULL,
                noFam        = NULL,
                rf           = NULL,
                LOD          = NULL,
                famInfo      = NULL,
                para         = NULL,
                LG           = NULL,
                LG_mat       = NULL,
                LG_mat_temp  = NULL,
                LG_pat       = NULL,
                LG_pat_temp  = NULL,
                ############################################
                ## function for mapping BI SNPs to maternal or paternal LGs
                mapBISnps = function(unmapped, parent, LODthres, nComp){
                  if(parent == "maternal")
                    newLGlist <- private$LG_mat
                  else if(parent == "paternal")
                    newLGlist <- private$LG_pat
                  nLG <- length(newLGlist)
                  ## Run algorithm for generating the linkage groups
                  noneMapped = FALSE
                  count = 0
                  while(!noneMapped){
                    noneMapped = TRUE
                    ## check that there are still SNPs remaining that need to be mapped
                    if(length(unmapped) == 0)
                      next
                    ## run the algorithm to map the SNPs
                    else{
                      for(snp in unmapped){
                        LODvalue = numeric(nLG)
                        for(lg in 1:nLG)
                          LODvalue[lg] <- mean(sort(private$LOD[snp,newLGlist[[lg]]],decreasing=T)[1:nComp],na.rm=T)
                        if(max(LODvalue) >= LODthres & sort(LODvalue, decreasing = T)[2] < LODthres){
                          count = count + 1
                          newLG <- which.max(LODvalue)
                          newLGlist[[newLG]] <- c(newLGlist[[newLG]], snp)
                          unmapped <- unmapped[-which(unmapped == snp)]
                          noneMapped = FALSE
                        }
                      }
                    }
                  }
                  return(newLGlist)
                }
              )
)

#### Make an unrelated population
### Make a full-sib family population
makePop.FS <- function(R6obj, pedfile, family=NULL, filter=list(MAF=0.05, MISS=0.2, BIN=0, DEPTH=5, PVALUE=0.01), inferSNPs = FALSE, perInfFam=1){
  
  ## Do some checks
  if(is.null(filter$MAF) || filter$MAF<0 || filter$MAF>1 || !is.numeric(filter$MAF)){
    warning("Minor allele filter has not be specified or is invalid. Setting to 5%:")
    filter$MAF <- 0.05
  }
  if(perInfFam <= 0.5 || perInfFam > 1)
    stop("The of the percentage of families which are informative for each SNP must greater than 50% and leass than equal to 100%.")
  if(is.null(filter$MISS) || filter$MISS<0 || filter$MISS>1 || !is.numeric(filter$MISS)){
    warning("Proportion of missing data filter has not be specified or is invalid. Setting to 20%:")
    filter$MISS <- 0.2
  }
  if(is.null(filter$BIN) || filter$BIN<0 || !is.finite(filter$BIN) || !is.numeric(filter$BIN)){
    warning("Minimum distance between adjacent SNPs is not specified or is invalid. Setting to 0:")
    filter$BIN <- 0 
  }
  if(is.null(filter$DEPTH) || filter$DEPTH<0 || is.infinite(filter$DEPTH) || !is.numeric(filter$DEPTH)){
    warning("Minimum depth on the parental genotypes filter has not be specified or is invalid. Setting to a depth of 5")
    filter$DEPTH <- 5
  }
  if(is.null(filter$PVALUE) || filter$PVALUE<0 || is.infinite(filter$PVALUE) || !is.numeric(filter$PVALUE)){
    warning("P-value for segregation test is not specified or invalid. Setting a P-value of 0.01:")
    filter$PVALUE <- 0.05
  }
  
  ## Define variables that will be used.
  cat("-------------\n")
  cat("Processing Data.\n\n")
  
  cat("Filtering criteria for removing SNPs:\n")
  cat("Minor allele frequency (MAF) < ", filter$MAF,"\n", sep="")
  cat("Percentage of missing genotypes > ", filter$MISS*100,"%\n", sep="")
  cat("Read depth associated with at least one parental genotype <= ", filter$DEPTH,"\n", sep="")
  cat("P-value for segregation test < ", filter$PVALUE,"%\n\n", sep="")
  ## Extract the private variables we want
  indID <- R6obj$.__enclos_env__$private$indID
  nSnps <- R6obj$.__enclos_env__$private$nSnps
  
  ## sort out the pedigree
  ped <- read.csv(pedfile, stringsAsFactors=F)
  famInfo = list()
  ## work out how many families there are
  parents <- unique(ped[c("Mother","Father")])
  # check that parents are in the file
  parents <- parents[which(apply(parents,1, function(x) all(x %in% ped$IndividualID))),]
  ## find the families to be used
  if(is.null(family)){
    family <- unique(ped$Family)[!is.na(unique(ped$Family))]
    family <- family[-which(nchar(family)==0)]
  }
  else{
    if(any(!(family %in% unique(ped$Family)[!is.na(unique(ped$Family))])))
      stop("Family missing from the pedigree file. Please check the family ID suppied or the pedigree file.")
  }
  ## Create each family
  for(fam in family){
    progIndx <- which(ped$Family == fam)
    famParents <- unique( ped[progIndx,c("Mother","Father")] )
    if(nrow(famParents) > 1)
      stop(paste0("Individuals with the same family name (Family ",fam,") have different different parents. Please check the pedigree file"))
    if(any(is.na(famParents)))
      stop(paste0("Family ",fam," has missing parents. Please check the pedigree file"))
    father <- as.integer(famParents["Father"])
    mother <- as.integer(famParents["Mother"])
    ## Create the pedigree
    famInfo[[as.character(fam)]]$progeny <- ped$SampleID[progIndx]
    famInfo[[as.character(fam)]]$parents <- list(Father=ped$SampleID[which(ped$IndividualID==father)],
                                                 Mother=ped$SampleID[which(ped$IndividualID==mother)])
    ## Check to see if there are any grandparents
    grandparents = ped[which(ped$IndividualID==father),c("Mother","Father")]
    if(nrow(grandparents) == 1){
      famInfo[[as.character(fam)]]$grandparents$paternalGrandFather <- ped$SampleID[which(ped$IndividualID == grandparents[,"Father"])]
      famInfo[[as.character(fam)]]$grandparents$paternalGrandMother <- ped$SampleID[which(ped$IndividualID == grandparents[,"Mother"])]
    }
    else if(nrow(grandparents) > 2)
      stop("Father of family ",fam," has mutiple parents. Please check the pedigree file.")
    grandparents = ped[which(ped$IndividualID==mother),c("Mother","Father")]
    if(nrow(grandparents) == 1){
      famInfo[[as.character(fam)]]$grandparents$maternalGrandFather <- ped$SampleID[which(ped$IndividualID == grandparents[,"Father"])]
      famInfo[[as.character(fam)]]$grandparents$maternalGrandMother <- ped$SampleID[which(ped$IndividualID == grandparents[,"Mother"])]
    }   
    else if(nrow(grandparents) > 2)
      stop("Mother of family ",fam," has mutiple parents. Please check the pedigree file.")
  }
  
  noFam <- length(famInfo)
  config_all <- config_infer_all <- nInd_all <- indx <- indID_all <- vector(mode = "list", length = noFam)
  genon_all <- ref_all <- alt_all <- vector(mode="list", length=noFam)
  
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
    ref <- R6obj$.__enclos_env__$private$ref[progIndx,]
    alt <- R6obj$.__enclos_env__$private$alt[progIndx,]
    
    ## Determine the segregation types of the loci
    genon_mum <- matrix(R6obj$.__enclos_env__$private$genon[mumIndx,], nrow=length(mumIndx), ncol=nSnps) 
    genon_dad <- matrix(R6obj$.__enclos_env__$private$genon[dadIndx,], nrow=length(mumIndx), ncol=nSnps)
    depth_mum <- matrix(R6obj$.__enclos_env__$private$ref[mumIndx,] +
                          R6obj$.__enclos_env__$private$alt[mumIndx,], nrow=length(mumIndx), ncol=nSnps)
    depth_dad <- matrix(R6obj$.__enclos_env__$private$ref[dadIndx,] +
                          R6obj$.__enclos_env__$private$alt[dadIndx,], nrow=length(mumIndx), ncol=nSnps)
    
    if(patgrandparents){
      genon_patgrandmum <- matrix(R6obj$.__enclos_env__$private$genon[patgrandmumIndx,], nrow=length(patgrandmumIndx), ncol=nSnps) 
      depth_patgrandmum <- matrix(R6obj$.__enclos_env__$private$ref[patgrandmumIndx,] +
                                    R6obj$.__enclos_env__$private$alt[patgrandmumIndx,], nrow=length(patgrandmumIndx), ncol=nSnps)
      genon_patgranddad <- matrix(R6obj$.__enclos_env__$private$genon[patgranddadIndx,], nrow=length(patgranddadIndx), ncol=nSnps) 
      depth_patgranddad <- matrix(R6obj$.__enclos_env__$private$ref[patgranddadIndx,] +
                                    R6obj$.__enclos_env__$private$alt[patgranddadIndx,], nrow=length(patgranddadIndx), ncol=nSnps)
    }
    if(matgrandparents){
      genon_matgrandmum <- matrix(R6obj$.__enclos_env__$private$genon[matgrandmumIndx,], nrow=length(matgrandmumIndx), ncol=nSnps) 
      depth_matgrandmum <- matrix(R6obj$.__enclos_env__$private$ref[matgrandmumIndx,] +
                                    R6obj$.__enclos_env__$private$alt[matgrandmumIndx,], nrow=length(matgrandmumIndx), ncol=nSnps)
      genon_matgranddad <- matrix(R6obj$.__enclos_env__$private$genon[matgranddadIndx,], nrow=length(matgranddadIndx), ncol=nSnps) 
      depth_matgranddad <- matrix(R6obj$.__enclos_env__$private$ref[matgranddadIndx,] +
                                    R6obj$.__enclos_env__$private$alt[matgranddadIndx,], nrow=length(matgranddadIndx), ncol=nSnps)
    }
    
    parHap_pat <- sapply(1:nSnps,function(x){
      x_p = genon_dad[,x]; d_p = depth_dad[,x]
      if(any(x_p==1,na.rm=T))
        return("AB")
      else if(sum(d_p) > filter$DEPTH){
        if(all(x_p==2, na.rm=T))
          return("AA")
        else if(all(x_p==0, na.rm=T))
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
        if(all(x_m==2, na.rm=T))
          return("AA")
        else if(all(x_m==0, na.rm=T))
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
    config[which(parHap_pat == "BB" & parHap_mat == "AB")] <- 5
    
    #### Segregation test to determine if the SNPs have been miss-classified
    seg_Dis <- sapply(1:nSnps,function(x){
      if(is.na(config[x]))
        return(NA)
      else{
        d = ref[,x] + alt[,x]
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
          d = ref[,x] + alt[,x]
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
    
    nInd_all[[fam]] <- nInd
    
    genon_all[[fam]] <- genon
    ref_all[[fam]] <- ref
    alt_all[[fam]] <- alt
  }
  
  if(noFam == 1){
    fam = 1
    
    ## Find all the SNPs to keep and subset the global variables
    indx_all <- indx[[fam]]
    
    #indx_all <- do.call("rbind",indx)
    #indx_all <- apply(indx_all, 2, any)
    
    genon_all[[fam]]     <- genon[,indx_all]
    ref_all[[fam]] <- ref[,indx_all]
    alt_all[[fam]] <- alt[,indx_all]
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
    tabInf_MI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(4,5) & indx[[y]]))))
    group.temp$MI <- setdiff(as.numeric(names(tabInf_MI)[which(tabInf_MI >= noInfoFam)]), group.temp$BI)
    # PI 
    tabInf_PI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(2,3) & indx[[y]]))))
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
      ref_all[[fam]] <- ref_all[[fam]][,indx_all]
      alt_all[[fam]] <- alt_all[[fam]][,indx_all]
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
    genon = genon_all, ref = ref_all, alt = alt_all, chrom = chrom_all, pos = pos_all,
    group = group, group_infer = group_infer, config = config_all, config_infer = config_infer_all,
    nInd = nInd_all, nSnps = sum(indx_all), noFam = noFam, indID = indID_all, SNP_Names = SNP_Names,
    masked=rep(FALSE,sum(indx_all)), famInfo=famInfo)
  )
  
  return(R6obj)
}


### Function for extending a vector to length n
extendVec <- function(vec, n){
  if(length(vec) == n)
    return(vec)
  else if (length(vec) < n){
    return(vec[1:n])
  }
  else{
    temp <- rep(NA,n)
    temp[1:length(vec)] <- vec
    return(temp)
  }
}
