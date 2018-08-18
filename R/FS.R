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
#' FS object
#' 
#' Class for storing RA data and associated functions for analysis with full-sib populations.
#' 
#' @usage
#' ## Create FS object
#' FSobj <- makeFS(RAobj, pedfile, family=NULL, 
#'                 filter=list(MAF=0.05, MISS=0.2, BIN=0, DEPTH=5, PVALUE=0.01))
#'
#' ## Functions (Methods) of FS object
#' FSobj$maskSNP(snps)
#' FSobj$rf_2pt(nClust=2)
#' FSobj$unmaskSNP(snps)
#' 
#' @details
#' An FS object is created from the \code{\link{makeFS}} function and contains RA data,
#' various statistics of the dataset that have been computed, and functions (methods)
#' for analyzing the data. Information in an FS object are specific to full-sib family populations.
#' 
#' @section Methods(Functions):
#' \describe{
#' \item{\code{\link{$maskSNP}}}{Mask SNPs from the data set.}
#' \item{\code{\link{$rf_2pt}}}{Compute the 2-point recombination fraction between all SNP pairs.}
#' \item{\code{\link{$unmaskSNP}}}{Unmask SNPs from the data set.}
#' }
#' @format NULL
#' @author Timothy P. Bilton
#' @name FS
#' @export

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
                  private$para      <- NULL
                },
                ### print methods
                 print = function(...){
                   if(private$noFam == 1){
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
                       cat("Linkage Group Summary:\n")
                       MI <- unlist(lapply(private$LG_mat, length))
                       PI <- unlist(lapply(private$LG_pat, length))
                       tab <- cbind(LG=1:(length(MI)+length(PI)),MI=c(MI,rep(0,length(PI))),PI=c(rep(0,length(MI)),PI))
                       prmatrix(tab, rowlab = rep("",nrow(tab)))
                     }
                     else{
                       cat("Linkage Group Summary:\n")
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
                       }
                       prmatrix(tab, rowlab = rep("",nrow(tab)))
                     }
                   }
                   else
                     cat("Not yet implemented")
                   return(self)
                 },
                #############################################################
                ## Function for removing SNPs from the linkage groups
                removeSNP = function(snps){
                  ## some checks
                  if( !is.vector(snps) || !is.numeric(snps)  || any(is.na(snps)) || !all(snps == round(snps)) || any(snps < 1) || any(snps > private$nSnps) )
                    stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                  snps <- unique(snps) ## make sure no double ups in SNPs
                  if(is.null(private$LG)){
                    if(is.null(private$LG_mat) || is.null(private$LG_pat))
                      stop("No linkage groups exist. Please use the $createLG function to create some linkage groups")
                    else{
                      ## remove SNP from the maternal linkage groups
                      for(lg in 1:length(private$LG_mat)){
                        if(any(private$LG_mat[[lg]] %in% snps))
                          private$LG_mat[[lg]] <- private$LG_mat[[lg]][-which(private$LG_mat[[lg]] %in% snps)]
                      }
                      ## remove SNP from the paternal linkage groups
                      for(lg in 1:length(private$LG_pat)){
                        if(any(private$LG_pat[[lg]] %in% snps))
                          private$LG_pat[[lg]] <- private$LG_pat[[lg]][-which(private$LG_pat[[lg]] %in% snps)]
                      }
                      ## check that there are no LGs that are now empty
                      empty_mat <- lapply(private$LG_mat, length) != 0
                      if(any(empty_mat))
                        private$LG_mat <- private$LG_mat[which(empty_mat)]
                      empty_pat <- lapply(private$LG_pat, length) != 0
                      if(any(empty_pat))
                        private$LG_pat <- private$LG_pat[which(empty_pat)]
                    }
                  }
                  else{
                    ## remove SNPs from the combined groups
                    for(lg in 1:length(private$LG))
                      private$LG[[lg]] <- private$LG[[lg]][-which(private$LG[[lg]] %in% snps)]
                    ## check for empty LGs
                    empty <- lapply(private$LG, length) != 0
                    if(any(empty))
                      private$LG <- private$LG[which(empty)]
                  }
                  return(invisible())
                },
                ## Function for removing linkage groups 
                removeLG = function(LG){
                  if(is.null(private$LG)){
                    if(is.null(private$LG_mat) || is.null(private$LG_pat))
                      stop("No linkage groups exist. Please use the $createLG function to create some linkage groups")
                    else{
                      nmat <- length(private$LG_mat)
                      npat <- length(private$LG_pat)
                      if(isValue(LG,min=1,max=nmat+npat))
                        stop(paste0("The LG indices must an integer vector between 1 and the number of linkage groups which is ",nmat+npat))
                      else{
                        private$LG_mat <- private$LG_mat[which(!(1:nmat %in% LG))]
                        private$LG_pat <- private$LG_pat[which(!((nmat+1:npat) %in% LG))]
                      }
                    }
                  }
                  else{
                    nLG <- length(private$LG)
                    if(isValue(LG,min=1,max=nLG))
                      stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nLG))
                    else{
                      private$LG <- private$LG[which(!(1:nLG %in% LG))]
                    }
                  }
                },
                ## function for merging linkage groups
                mergeLG = function(LG, mergeTo=NULL){
                  ## checks
                  if(!is.null(mergeTo) && (!is.character(mergeTo) || length(mergeTo) != 1 || !(mergeTo %in% c("maternal","paternal","none"))))
                    stop("Argument for which parent group is to be merged to is invalid")
                  if(is.null(private$LG)){
                    ## Check that we have linkage groups
                    if(is.null(private$LG_mat) || is.null(private$LG_pat))
                      stop("No linkage groups exist. Please use the $createLG function to create some linkage groups")
                    else{
                      nmat <- length(private$LG_mat)
                      npat <- length(private$LG_pat)
                      ## check input
                      if(isValue(LG,min=1,max=nmat+npat))
                        stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nmat+npat))
                      LG <- unique(LG)
                      if(length(LG) == 1)
                        stop("Only one unique linkage group supplied. Need at least two linkage groups to merge")
                      matgroups <- LG[which(LG %in% 1:nmat)]
                      patgroups <- LG[which(LG %in% (nmat + 1:npat))] - nmat
                      ## merge maternals 
                      if(length(matgroups) > 1)
                        private$LG_mat[[matgroups[1]]] <- unlist(private$LG_mat[matgroups])
                      ## merge paternals
                      if(length(patgroups) > 1)
                        private$LG_pat[[patgroups[1]]] <- unlist(private$LG_pat[patgroups])
                      ## merge between parent meiosis
                      if(length(matgroups) > 0 & length(patgroups) > 0 & (is.null(mergeTo) || mergeTo != "none")){
                        if(is.null(mergeTo)){
                          mergeTo <- menu(c("maternal","paternal"),title="Merging maternal and paternal LGs...\nWhich parental meiosis to merge groups?")
                          if(mergeTo == 0)
                            stop("Linkage groups between parental meiosis were NOT merged.")
                          mergeTo <- switch(mergeTo, "maternal", "paternal")
                        }
                        ## merge groups
                        mergedLG <- unlist(c(private$LG_mat[matgroups[1]],private$LG_pat[patgroups[1]]))
                        if(mergeTo == "maternal"){
                          private$LG_mat[[matgroups[1]]] <- mergedLG
                          private$LG_mat[matgroups[-1]] <- NULL
                          private$LG_pat[patgroups] <- NULL
                        } else if(mergeTo == "paternal"){
                          private$LG_pat[[patgroups[1]]] <- mergedLG
                          private$LG_mat[matgroups] <- NULL
                          private$LG_pat[patgroups[-1]] <- NULL
                        }
                      } else{
                        ## remove the LGs that were merged
                        if(length(matgroups) > 1)
                          private$LG_mat[matgroups[-1]] <- NULL
                        if(length(patgroups) > 1)
                          private$LG_pat[patgroups[-1]] <- NULL
                      }
                    }
                  } else{
                    nLG <- length(private$LG)
                    ## check input
                    if(isValue(LG,min=1,max=nLG))
                      stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nLG))
                    LG <- unique(LG)
                    if(length(LG) == 1)
                      stop("Only one unique linkage group supplied. Need at least two linkage groups to merge")
                    ## merge LGs into the first
                    private$LG[[LG[1]]] <- unlist(private$LG[LG])
                    ## drop remaining groups
                    private$LG[LG[-1]] <- NULL
                  }
                  return(invisible())
                },
                ## function for masking SNPs
                maskSNP = function(snps){
                  ## Inout checks
                  if( !is.vector(snps) || !is.numeric(snps) || any(is.na(snps)) || !all(snps == round(snps)) || any(snps < 1) || any(snps > private$nSnps) )
                    stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                  ## mask SNPs
                  private$masked[snps] <- TRUE
                  return(invisible(self))
                },
                ## function for unmasking SNPs
                unmaskSNP = function(snps){
                  if( !is.vector(snps) || !is.numeric(snps) || any(is.na(snps)) || !all(snps == round(snps)) || any(snps < 1) || any(snps > private$nSnps) )
                    stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                  private$masked[snps] <- FALSE
                  return(self)
                },
                #####################################################################
                ## Function for computing the 2-point rf estimates
                rf_2pt = function(nClust=4, err=FALSE){
                  ## do some checks
                  if(!is.numeric(nClust) || !is.vector(nClust) || length(nClust)!=1 || nClust < 0 || is.infinite(nClust) )
                    stop("Number of clusters for parallelization needs to be a positive finite integer number")
                  if(!is.vector(err) || !is.logical(err) || length(err) != 1)
                    stop("Argument specifying whether error parameters are to be estimated is invlaid.")
                  ## If there is only one family
                  if(private$noFam == 1)
                    mat <- rf_2pt_single(private$ref[[1]], private$alt[[1]],
                                         private$config[[1]], private$config_infer[[1]],
                                         private$group, private$group_infer,
                                         nClust, private$nInd, err=err)
                  else
                    stop("Multiple families have yet to be implemented")
                    #mat <- rf_2pt_multi(private$ref, private$alt,
                    #                    private$config,private$group, nClust, private$noFam)
                  ## Save the results to the object
                  private$rf <- mat$rf
                  private$LOD <- mat$LOD
                  return(invisible())
                },
                ## Function for creating linkage groups
                createLG = function(parent="both", LODthres=10, nComp=10){
                  ## Do some checks
                  if(is.null(private$rf) || is.null(private$LOD))
                    stop("Recombination fractions and LOD scores have not been computed.\nUse $rf_2pt() to compute the recombination fractions and LOD scores.")
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop("parent argument is not a string of length one or is invalid:
Please select one of the following:
   maternal: Only MI SNPs
   paternal: Only PI SNPs
   both:     MI and PI SNPs")
                  if(!is.numeric(LODthres) || !is.vector(LODthres) || length(LODthres) != 1 || is.na(nComp) ||  LODthres < 0 || !is.finite(LODthres))
                    stop("The LOD threshold (argument 2) needs to be a finite numeric number.")
                  if(!is.numeric(nComp) || !is.vector(nComp) || length(nComp) != 1 || is.na(nComp) || nComp < 0 || !is.finite(nComp) ||
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
                addBIsnps = function(LODthres=10, nComp=10){

                  ## Do some checks
                  if(is.null(private$rf) || is.null(private$LOD))
                    stop("Recombination fractions and LOD scores have not been computed.\nUse rf_2pt() to compute the recombination fractions and LOD scores.")
                  if(is.null(private$LG_pat) && is.null(private$LG_mat))
                    stop("There are no existing linkage groups. Use $createLG() method to create linkage groups.")
                  ##if(parent == "maternal" && length(private$LG_mat) == 0)
                  #  stop("There are no maternal linkage groups. Please use 'createLG' function to create the linkage groups first")
                  #else if(parent == "paternal" && length(private$LG_pat) == 0)
                  #  stop("There are no paternal linkage groups. Please use 'createLG' function to create the linkage groups first")
                  if(!is.numeric(LODthres) || !is.vector(LODthres) || length(LODthres) != 1 || LODthres < 0 || !is.finite(LODthres))
                    stop("The LOD threshold (argument 1) needs to be a finite positive numeric number.")
                  if(!is.numeric(nComp) || !is.vector(nComp) || length(nComp) != 1 || nComp < 0 || !is.finite(nComp) ||
                     round(nComp) != nComp)
                    stop("The number of comparsion points (argument 2) needs to be a finite positive integer number.")

                  ## Find the unmapped loci
                  unmapped <- sort(unlist(private$group$BI,private$group_infer$BI))
                  ## Remove masked SNPs
                  unmapped <- unmapped[which(!private$masked[unmapped])]
                  ## check that there are SNPs to map
                  if(length(unmapped) == 0)
                    stop("There are no SNPs remaining that are unmapped")
                  
                  ## if maternal groups, map the BIs
                  if(!is.null(private$LG_mat))
                    newLGlist_mat <- private$mapBISnps(unmapped, parent="maternal", LODthres=LODthres, nComp=nComp)
                  ## if paternal groups, map the BIs
                  if(!is.null(private$LG_pat))
                    newLGlist_pat <- private$mapBISnps(unmapped, parent="paternal", LODthres=LODthres, nComp=nComp)
                  
                  if(exists("newLGlist_mat") & exists("newLGlist_pat")){
                    
                    ## work out which LGs go together
                    indx_mat <- sapply(unmapped, function(x){ which(unlist(lapply(newLGlist_mat, function(y) any(x == y))))}, simplify=F)
                    indx_pat <- sapply(unmapped, function(x){ which(unlist(lapply(newLGlist_pat, function(y) any(x == y))))}, simplify=F)
                    mappedBI <- which(sapply(1:length(unmapped), function(y) length(indx_mat[[y]]) != 0 &  length(indx_pat[[y]])))
                    nmat <- length(private$LG_mat)
                    npat <- length(private$LG_pat)
                    tabBI <- table(unlist(factor(indx_mat[mappedBI], levels=1:nmat)),
                                   unlist(factor(indx_pat[mappedBI], levels=1:npat)))
                    finish <- FALSE
                    mappedLG_mat <- mappedLG_pat <- numeric(0)
                    while(!finish){
                      if(length(mappedLG_mat) == nmat || length(mappedLG_pat) == npat){
                        finish = TRUE
                        next
                      }
                      maxv <- max(tabBI[which(!(1:nmat %in% mappedLG_mat)),which(!(1:npat %in% mappedLG_pat))], na.rm=T)
                      if(maxv == 0){
                        finish = TRUE
                        next
                      }
                      temp <- which(tabBI == maxv, arr.ind=T)[1,]
                      tabBI[t(temp)] <- NA
                      mappedLG_mat <- c(mappedLG_mat,(1:nmat)[temp[1]])
                      mappedLG_pat <- c(mappedLG_pat,(1:npat)[temp[2]])
                    }
                    nmapped <- length(mappedLG_mat)
                    ## Check whether there are still any lingering LGs that could be mapped
                    if(unique(length(mappedLG_mat)) < nmat){
                      tomap <- which(!(1:nmat %in% mappedLG_mat))
                      if(max(tabBI[tomap,], na.rm=T) > 0){
                        for(matlg in tomap){
                          patlg <- which.max(tabBI[matlg,])
                          mappedLG_mat <- c(mappedLG_mat,matlg)
                          mappedLG_pat <- c(mappedLG_pat,patlg)
                          tabBI[matlg,patlg] <- NA
                        }
                      }
                    }
                    if(unique(length(mappedLG_pat)) < npat){
                      tomap <- which(!(1:npat %in% mappedLG_pat))
                      if(max(tabBI[,tomap], na.rm=T) > 0){
                        for(patlg in tomap){
                          matlg <- which.max(tabBI[,patlg])
                          mappedLG_mat <- c(mappedLG_mat,matlg)
                          mappedLG_pat <- c(mappedLG_pat,patlg)
                          tabBI[matlg,patlg]
                        }
                      }
                    }
                    ## Now combine groups with only the BIs mapped to same groups
                    LGmerged <- list()
                    for(lg in 1:nmapped){
                      LGmerged[[lg]] <- c(newLGlist_mat[[mappedLG_mat[lg]]], newLGlist_pat[[mappedLG_pat[lg]]])
                    }
                    ## check for extra groups that have been mapped
                    if(nmapped != length(mappedLG_mat)){
                      for(lg in (nmapped+1):length(mappedLG_mat)){
                        matlg <- which(mappedLG_mat[1:nmapped] == mappedLG_mat[lg])
                        patlg <- which(mappedLG_pat[1:nmapped] == mappedLG_pat[lg])
                        if(length(matlg) == 0)
                          LGmerged[[patlg]] <- c(LGmerged[[patlg]],newLGlist_mat[[mappedLG_mat[lg]]])
                        else if(length(patlg) == 0)
                          LGmerged[[matlg]] <- c(LGmerged[[matlg]],newLGlist_mat[[mappedLG_pat[lg]]])
                      }
                    }
                    ## Clean-up: remove BI SNPs which only mapped to one LG and remove duplicates
                    for(lg in 1:length(LGmerged)){
                      bi_indx <- LGmerged[[lg]] %in% private$group$BI
                      bisnps <- LGmerged[[lg]][which(bi_indx)]
                      LGmerged[[lg]] <- c(LGmerged[[lg]][which(!bi_indx)],unique(bisnps[duplicated(bisnps)]))
                    }
                    private$LG <- LGmerged
                  }
                  else if(exists("newLGlist_mat")){
                    private$LG <- newLGlist_mat
                  }
                  else if(exists("newLGlist_pat")){
                    private$LG <- newLGlist_pat
                  }
                  return(invisible())
                },
                ## Function for ordering linkage groups
                orderLG = function(chrom = NULL, mapfun = "morgan", weight="LOD2", ndim=2, spar=NULL){
                  ## do some checks
                  if(is.null(private$LG))
                    stop("There are no combined linkage groups. Please use the $addBISNPs to create combined LGs")
                  if(is.null(chrom))
                    chrom <- 1:length(private$LG)
                  else if(!is.vector(chrom) && !is.numeric(chrom) && any(chrom < 0) && chrom > length(private$LG))
                    stop("Invalid chromosome input")
                  if(!(mapfun %in% c("morgan","haldane","kosambi")))
                    stop("Unknown mapping function")
                  if(!(weight %in% c("LOD","LOD2","none")))
                     stop("Unknown weighting function")
                  nChr = length(chrom)
                  LG <- vector(mode='list', length=nChr)
                  ## order the linkage groups
                  for(chr in chrom){
                    ind <- private$LG[[chr]]
                    ind_mi <- which(ind %in% private$group$MI)
                    ind_pi <- which(ind %in% private$group$PI)
                    ## set up the weighting matrix
                    if(weight == "LOD")
                      wmat <- matrix(private$LOD[ind,ind],nrow=length(ind), ncol=length(ind))
                    else if (weight == "LOD2")
                      wmat <- matrix((private$LOD[ind,ind])^2,nrow=length(ind), ncol=length(ind))
                    else if (weight == "none")
                      wmat <- matrix(1,nrow=length(ind), ncol=length(ind))
                    ## set the weights of the PI and MI combinations to zero.
                    wmat[ind_mi, ind_pi] <- 0
                    wmat[ind_pi, ind_mi] <- 0
                    ## set of the distance matrix
                    dmat <- mfun(private$rf[ind,ind], fun = mapfun, centiM=FALSE)
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
                    private$LG[[chr]] <- ind[pcurve$ord]
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
                plotLG = function(mat=c("rf"), parent, LG=NULL, filename=NULL, names=NULL, chrS=2, lmai=2, chrom=T){
                  ## do some checks
                  if(!is.vector(mat) || !is.character(mat) || length(mat) != 1 || !(mat %in% c('rf','LOD')))
                    stop("Argument specifying which matrix to plot (argument 1) must be either 'rf' or 'LOD'")
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop("parent argument is not a string of length one or is invalid:
Please select one of the following:
  maternal: Add BI SNPs to MI LGs
  paternal: Add BI SNPs to PI LGs
  both:     Add BI SNPs to both MI and PI LGs")
                  
                  if(private$noFam == 1){
                    ## Work out which LGs list to use
                    if(is.null(private$LG)){
                      if(!is.null(private$LG_pat) || !is.null(private$LG_mat))
                        LGlist <- c(private$LG_mat,private$LG_pat)
                      else
                        stop("No linkage groups are available to be plotted. Please use the $createLG function to create some linkage groups")
                    }
                    else
                      LGlist <- private$LG
                    
                    ## Work out if we want a subset of the LGs
                    if(!is.null(LG)){
                      LGlist <- LGlist[LG]
                      names(LGlist) <- LG
                    }
                    else
                      names(LGlist) <- 1:length(LGlist)
                    
                    
                    ## Check which type of SNPs we are plotting
                    if(parent == "maternal"){
                      mi_ind <- lapply(LGlist, function(x) x[which(x %in% c(private$group$MI, private$group$BI))])
                      LGlist <- mi_ind[which(unlist(lapply(mi_ind, length))!=0)]
                    }
                    else if (parent == "paternal"){
                      pi_ind <- lapply(LGlist, function(x) x[which(x %in% c(private$group$PI, private$group$BI))])
                      LGlist <- pi_ind[which(unlist(lapply(pi_ind, length))!=0)]
                    }
  
                    ## Sort out the matrix
                    if(mat == "rf")
                      temprf <- private$rf
                    else if(mat == "LOD")
                      temprf <- private$LOD
                    b <- ncol(temprf) + 1
                    temprf <- cbind(temprf,rep(NA,b-1))
                    temprf <- rbind(temprf,rep(NA,b))
                    chrom.ind <- unlist(lapply(LGlist, function(x) c(x,b)), use.names = FALSE)
                    chrom.ind <- chrom.ind[-length(chrom.ind)]
                    temprf <- temprf[chrom.ind, chrom.ind]
                    
                    ## Plot the matrix
                    nn <- length(chrom.ind)
                    b_indx <- chrom.ind == b
                    chrom.ind[which(!b_indx)] <- paste0(chrom.ind[which(!b_indx)]," (", rep(names(LGlist), lapply(LGlist,length)),")") 
                    chrom.ind[which(b_indx)] <- rep("Break",length(LGlist)-1)
                    hovertext <- matrix(paste(matrix(paste0("x: ",chrom.ind), nrow=nn, ncol=nn), 
                          matrix(paste0("x: ",chrom.ind), nrow=nn, ncol=nn, byrow=T), paste0("rf: ",round(temprf,4)), sep="<br>"),
                          nrow=nn, ncol=nn)
                    ax <- list(visible=FALSE)
                    # suppress warnings  
                    storeWarn<- getOption("warn")
                    options(warn = -1) 
                    ## produce the plotly plot
                    p <- plotly::plot_ly(z=temprf, type="heatmap", showscale=F, hoverinfo="text",
                            text=hovertext, colors=heat.colors(100)) %>% 
                      plotly::add_segments(x=which(b_indx)-1,xend=which(b_indx)-1,y=0,yend=nn, line=list(color="black"),  showlegend=F) %>%
                      plotly::add_segments(y=which(b_indx)-1,yend=which(b_indx)-1,x=0,xend=nn, line=list(color="black"),  showlegend=F) %>%
                      plotly::layout(margin=list(l=0,r=0,t=0,b=0), xaxis=ax, yaxis=ax)
                    print(p)
                    options(warn = storeWarn) 
                  }
                  else
                    stop("Not yet implemented.")
                  return(invisible())
                },
                ## Function for plotting chromosome in their original ordering
                plotChr = function(mat=c("rf"), parent = "maternal", filename=NULL, chrS=2, lmai=2){
                  ## do some checks
                  if(!is.vector(mat) || !is.character(mat) || length(mat) != 1 || !(mat %in% c('rf','LOD')))
                    stop("Argument specifying which matrix to plot (argument 1) must be either 'rf' or 'LOD'")
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop("parent argument is not a string of length one or is invalid:
Please select one of the following:
  maternal: Add BI SNPs to MI LGs
  paternal: Add BI SNPs to PI LGs
  both:     Add BI SNPs to both MI and PI LGs")
                  
                  ## workout which SNPs to plot to plot
                  if(private$noFam == 1){
                    ## workout the indices for the chromosomes
                    names <- unique(private$chrom)
                    if(parent == "maternal")
                      LG <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,4,5))), simplify=F)
                    else if(parent == "paternal")
                      LG <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,2,3))), simplify=F)
                    else if(parent == "both")
                      LG <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,2,3,4,5))), simplify=F)
                    ## plot the chromsomes rf info
                    if(mat == "rf")
                      plotLG(mat=private$rf, LG=LG, filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=T)
                    else if(mat == "LOD")
                      plotLG(mat=private$LOD, LG=LG, filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=T)
                    else
                      stop("Matrix to be plotted not found.") ## shouldn't get here
                    return(invisible())
                  }
                  else{
                    stop("not implemented yet")
                  }
                },
                ## compare order with original and ordered markers
                plotSyn = function(){
                  LGorder <- unlist(private$LG)
                  orgOrder <- sort(LGorder)
                  LGbreaks <- cumsum(unlist(lapply(private$LG, length))) + 0.5
                  LGbreaks <- LGbreaks[-length(LGbreaks)]
                  chrBreaks <- which(diff(private$chrom)==1) + 0.5
                  
                  plot(orgOrder,LGorder, pch=20,cex=0.8, xaxt="n", yaxt="n",ylab="Assembly Ordering", xlab="Linkage Group Ordering", 
                       ylim=c(min(orgOrder),max(orgOrder)), xlim=c(min(LGorder), max(LGorder)))
                  abline(h=chrBreaks)
                  abline(v=LGbreaks)
                  mtext(text = unique(private$chrom[orgOrder]), side = 2,
                        at = apply(cbind(c(min(orgOrder),chrBreaks),c(chrBreaks,max(orgOrder))),1,mean))
                  mtext(text = 1:length(private$LG), side = 1,
                        at = apply(cbind(c(min(LGorder),LGbreaks),c(LGbreaks,max(LGorder))),1,mean))
                }, 
                ## Function for computing the rf's for each chromosome 
                rf_est = function(chr=NULL, init_r=0.01, ep=0.001, method="optim", sexSpec=F, seqErr=T, mapped=T){
                  ## do some checks
                  if( !is.null(init_r) & !is.numeric(init_r) )
                    stop("Starting values for the recombination fraction needs to be a numeric vector or integer or a NULL object")
                  if( (length(ep) != 1 || !is.numeric(ep) || (ep <= 0 | ep >= 1)) )
                    stop("Value for the error parameters needs to be a single numeric value in the interval (0,1) or a NULL object")
                  ## for existing chromosome orders
                  nChr <- length(private$LG)
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
                    private$para <- list(OPGP=vector(mode = "list",length = nChr),rf_p=vector(mode = "list",length = nChr),
                                         rf_m=vector(mode = "list",length = nChr),ep=vector(mode = "list",length = nChr),
                                         loglik=vector(mode = "list",length = nChr))
                    for(i in chr){
                      cat("Chromosome: ",i,"\n")
                      indx_chr <- which((private$chrom == i) & !private$masked)
                      ref_temp <- lapply(private$ref, function(x) x[,indx_chr])
                      alt_temp <- lapply(private$alt, function(x) x[,indx_chr])
                      ## estimate OPGP's
                      if(!any(names(private$para$OPGP) != i) || sapply(1:private$noFam, function(x)
                        !is.null(private$para$OPGP[i][[x]]) && (length(private$para$OPGP[i][[x]]) != ncol(ref_temp[[x]])))){
                        tempOPGP <- list()
                        for(fam in 1:private$noFam){
                          tempOPGP <- c(tempOPGP,list(as.integer(infer_OPGP_FS(ref_temp[[fam]],alt_temp[[fam]],private$config[[fam]][indx_chr], method="EM"))))                        
                          }
                        private$para$OPGP <- tempOPGP
                      }
                      ## estimate the rf's
                      MLE <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP,
                                         sexSpec=sexSpec, seqErr=seqErr, method=method)
                      if(sexSpec){
                        private$para$rf_p[i]   <- list(MLE$rf_p)
                        private$para$rf_m[i]   <- list(MLE$rf_m)
                      } else{
                        private$para$rf_p[i]   <- list(MLE$rf)
                        private$para$rf_m[i]   <- list(MLE$rf)
                      }
                      private$para$ep[i]    <- list(list(MLE$ep))
                      private$para$loglik[i]<- list(MLE$loglik)
                    }
                  }
                  else{
                    ## Check the input
                    if((length(private$LG) == 0) && length(is.null(private$LG) == 0))
                      stop("No linkage groups have been formed.")
                    if(is.null(chr)) # compute distances for all chromosomes
                      chr <- 1:nChr
                    else if(!is.vector(chr) || chr < 0 || chr > length(private$LG) || round(chr) != chr)
                      stop("Chromosomes input must be an integer vector between 1 and the number of linkage groups")
                    ## Compute rf's
                    cat("Computing recombination fractions:\n")
                    private$para <- list(OPGP=vector(mode = "list",length = nChr),rf_p=vector(mode = "list",length = nChr),
                                         rf_m=vector(mode = "list",length = nChr),ep=vector(mode = "list",length = nChr),
                                         loglik=vector(mode = "list",length = nChr))
                    for(i in chr){
                      cat("Chromosome: ",i,"\n")
                      indx_chr <- private$LG[[chr]]
                      ref_temp <- lapply(private$ref, function(x) x[,indx_chr])
                      alt_temp <- lapply(private$alt, function(x) x[,indx_chr])
                      ## estimate OPGP's
                      if(!any(names(private$para$OPGP) != i) || sapply(1:private$noFam, function(x)
                        !is.null(private$para$OPGP[i][[x]]) && (length(private$para$OPGP[i][[x]]) != ncol(ref_temp[[x]])))){
                        tempOPGP <- list()
                        for(fam in 1:private$noFam){
                          tempOPGP <- c(tempOPGP,list(as.integer(infer_OPGP_FS(ref_temp[[fam]],alt_temp[[fam]],private$config[[fam]][indx_chr], method="EM"))))                        
                        }
                        private$para$OPGP[i] <- tempOPGP
                      }
                      ## estimate the rf's
                      MLE <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                       sexSpec=sexSpec, seqErr=seqErr, method=method)
                      if(sexSpec){
                        private$para$rf_p[i]   <- list(MLE$rf_p)
                        private$para$rf_m[i]   <- list(MLE$rf_m)
                      } else{
                        private$para$rf_p[i]   <- list(MLE$rf)
                        private$para$rf_m[i]   <- list(MLE$rf)
                      }
                      private$para$ep[i]       <- list(MLE$ep)
                      private$para$loglik[i]   <- list(MLE$loglik)
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
