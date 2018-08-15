#' @export FS

### R6 class for creating a data format for a full-sib families
FS <- R6Class("FS",
              inherit = RA,
              public = list(
                ## variables
                LG_mat    = NULL,
                LG_pat    = NULL,
                seg_group = NULL,
                seg_num   = NULL,
                para      = NULL,
                ## initialize function
                initialize = function(R6obj){
                  self$LG_mat       <- list()
                  self$LG_pat       <- list()
                  self$seg_group    <- list()
                  self$seg_num      <- list()
                  self$para         <- list(OPGP=list(),rf_p=list(),rf_m=list(),ep=list(),loglik=list())
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
                  private$masked    <- rep(FALSE, R6obj$.__enclos_env__$private$nSnps)
                  private$famInfo   <- R6obj$.__enclos_env__$private$famInfo
                },
                #############################################################
                ## Function for removing SNPs from the linkage groups
                removeSNP = function(indx){
                  ## some checks
                  if( !is.vector(indx) || !is.numeric(indx) || !all(indx == round(indx)) || any(indx < 1) || any(indx > private$nSnps) )
                    stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                  ## remove SNP from the maternal linkage groups
                  for(lg in 1:length(self$LG_mat)){
                    if(any(self$LG_mat[[lg]] %in% indx))
                      self$LG_mat[[lg]] <- self$LG_mat[[lg]][-which(self$LG_mat[[lg]] %in% indx)]
                  }
                  ## remove SNP from the paternal linkage groups
                  for(lg in 1:length(self$LG_pat)){
                    if(any(self$LG_pat[[lg]] %in% indx))
                      self$LG_pat[[lg]] <- self$LG_pat[[lg]][-which(self$LG_pat[[lg]] %in% indx)]
                  }
                  return(invisible())
                },
                ## function for masking SNPs
                maskSNP = function(indx){
                  if( !is.vector(indx) || !is.numeric(indx) || !all(indx == round(indx)) || any(indx < 1) || any(indx > private$nSnps) )
                    stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                  private$masked[indx] <- TRUE
                  return(invisible())
                },
                ## function for unmasking SNPs
                unMaskSNP = function(indx){
                  if( !is.vector(indx) || !is.numeric(indx) || !all(indx == round(indx)) || any(indx < 1) || any(indx > private$nSnps) )
                    stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                  private$masked[indx] <- FALSE
                  return(invisible())
                },
                
                ## Function for computing the 2-point rf estimates
                rf_2pt = function(nClust=4){
                  ## do some checks
                  if(!is.numeric(nClust) || length(nClust)!=1 || nClust < 0 || is.infinite(nClust) )
                    stop("Number clusters specifies for the parallelization needs to be a positive finite numeric number")
                  
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
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal")))
                    stop("parent argument is not a string of of length one or is incorrect:
                         Please select one of the following:
                         maternal: Only MI SNPs
                         paternal: Only PI SNPs")
                  if(!is.numeric(LODthres) || LODthres < 0 || !is.finite(LODthres) || length(LODthres) != 1)
                    stop("The LOD threshold (argument 2) needs to be a finite numeric number.")
                  if(!is.integer(nComp) || nComp < 0 || !is.finite(nComp) || length(nComp) != 1)
                    stop("The number of comparsion points (argument 3) needs to be a finite numeric number.")
                  
                  ## Create the groups
                  if(parent == "maternal")
                    self$LG_mat <- createLG(private$group, private$LOD, "maternal", LODthres, nComp, masked=private$masked)
                  else if(parent == "paternal")
                    self$LG_pat <- createLG(private$group, private$LOD, "paternal", LODthres, nComp, masked=private$masked)
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
                plotChr = function(mat=c("rf"), parent = "maternal", filename=NULL, chrS=2, lmai=2){
                  if(private$noFam == 1){
                    ## workout the indices for the chromosomes
                    names <- unique(private$chrom)
                    if(parent == "maternal")
                      LG <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,4,5))), simplify=F)
                    if(parent == "paternal")
                      LG <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,2,3))), simplify=F)
                    ## plot the chromsomes rf infom
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
                rf_est = function(chr=NULL, init_r=0.01, ep=0.001, method="optim", sexSpec=F, seqErr=T, nThreads=0){
                  ## do some checks
                  if( !is.null(init_r) & !is.numeric(init_r) )
                    stop("Starting values for the recombination fraction needs to be a numeric vector or integer or a NULL object")
                  if( (length(ep) != 1 || !is.numeric(ep) || (ep <= 0 | ep >= 1)) )
                    stop("Value for the error parameters needs to be a single numeric value in the interval (0,1) or a NULL object")
                  ## for existing chromosome orders
                  if((length(self$LG_mat) == 0) && length(is.null(self$LG_pat)==0)){
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
                          tempOPGP <- c(tempOPGP,list(as.integer(infer_OPGP_FS(ref_temp[[fam]],alt_temp[[fam]],private$config[[fam]][indx_chr], method="EM", nThreads=nThreads))))
                          }
                        self$para$OPGP[i] <- tempOPGP
                      }
                      ## estimate the rf's
                      cat(">>> Calling rf_est_FS", nThreads, "\n")
                      MLE <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=self$para$OPGP[i],
                                         sexSpec=sexSpec, seqErr=seqErr, method=method, nThreads=nThreads)
                      if(sexSpec){
                        self$para$rf_p[i]   <- MLE$rf_p
                        self$para$rf_m[i]   <- MLE$rf_m
                      } else{
                        self$para$rf_p[i]   <- list(MLE$rf)
                        self$para$rf_m[i]   <- list(MLE$rf)
                      }
                      self$para$ep[i]     <- MLE$ep
                      self$para$loglik[i] <- MLE$loglik
                    }
                  }
                  else{
                    stop("Yet to be implemented with grouping present")
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
                ## Function for updating the private variables
                updatePrivate = function(List){
                  if(!is.null(List$genon))
                    private$genon        = List$genon
                  if(!is.null(List$ref))
                    private$ref          = List$ref
                  if(!is.null(List$alt))
                    private$alt          = List$alt
                  if(!is.null(List$chrom))
                    private$chrom        = List$chrom
                  if(!is.null(List$pos))
                    private$pos          = List$pos
                  if(!is.null(List$SNP_Names))
                    private$SNP_Names    = List$SNP_Names
                  if(!is.null(List$indID))
                    private$indID        = List$indID
                  if(!is.null(List$nSnps))
                    private$nSnps        = List$nSnps
                  if(!is.null(List$nInd))
                    private$nInd         = List$nInd
                  if(!is.null(List$config))
                    private$config       = List$config
                  if(!is.null(List$config_infer))
                    private$config_infer = List$config_infer
                  if(!is.null(List$group))
                    private$group        = List$group
                  if(!is.null(List$group_infer))
                    private$group_infer  = List$group_infer
                  if(!is.null(List$noFam))
                    private$noFam        = List$noFam
                  if(!is.null(List$rf))
                    private$rf           = List$rf
                  if(!is.null(List$LOD))
                    private$LOD          = List$LOD
                  if(!is.null(List$masked))
                    private$masked       = List$masked
                  if(!is.null(List$famInfo))
                    private$famInfo      = List$famInfo
                }
              )
              )
