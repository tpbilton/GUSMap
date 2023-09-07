##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2020 Timothy P. Bilton
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
#' BC object
#' 
#' Class for storing BC data and associated functions for analysis of backcross family populations.
#' 
#' @usage
#' ## Create BC object
#' BCobj <- makeBC(RAobj, pedfile, family=NULL, inferSNPs=FALSE,
#'                 filter=list(MAF=0.05, MISS=0.2, BIN=100, DEPTH=5, PVALUE=0.01, MAXDEPTH=500))
#'
#' ## Functions (Methods) of an BC object
#' BCobj$addBIsnps(LODthres = 10, nComp = 10)
#' BCobj$iaddSNPs(LODthres = 10, nComp = 10)
#' BCobj$computeMap(chrom = NULL, init_r = 0.01, ep = 0.001, method = NULL, err = T, mapped = T, nThreads = 1, rfthres = 0.1)
#' BCobj$createLG(parent = "both", LODthres = 10, nComp = 10, reset = FALSE)
#' BCobj$maskSNP(snps)
#' BCobj$mergeLG(LG, where = NULL, mergeTo = NULL)
#' BCobj$orderLG(chrom = NULL, mapfun = "haldane", weight = "LOD2", ndim = 30, spar = NULL)
#' BCobj$plotChr(parent = "maternal", mat = "rf", filename = NULL, chrS = 2, lmai = 0.25)
#' BCobj$plotLG(parent  = "both", LG = NULL, mat = "rf", filename = NULL, interactive = FALSE, what = NULL)
#' BCobj$plotLM(LG = NULL, fun = "haldane", col = "black")
#' BCobj$plotSyn()
#' BCobj$print(what = NULL, ...)
#' BCobj$removeLG(LG, where = NULL)
#' BCobj$removeSNP(snps, where = NULL)
#' BCobj$rf_2pt(nClust = 2, err = FALSE)
#' BCobj$unmaskSNP(snps)
#' BCobj$writeLM(file, direct = "./", LG = NULL, what = NULL, inferGeno = TRUE)
#' 
#' @details
#' An BC object is created from the \code{\link{makeBC}} function and contains RA data,
#' various statistics of the dataset that have been computed, and functions (methods)
#' for analyzing the data. Information in an BC object are specific to backcross family populations.
#' 
#' @section Methods(Functions):
#' \describe{
#' \item{\code{\link{$addBIsnps}}}{Add Both-Informative (BI) snps to the maternal and/or
#' paternal linkage groups.}
#' \item{\code{\link{$addSNPs}}}{Add MI, PI and/or SI SNPs to the existing linkage groups}
#' \item{\code{\link{$createLG}}}{Create linkage group(s).}
#' \item{\code{\link{$computeMap}}}{Compute linkage maps for a given marker order.}
#' \item{\code{\link{$maskSNP}}}{Mask SNP(s) in the dataset.}
#' \item{\code{\link{$mergeLG}}}{Merge linkage groups.}
#' \item{\code{\link{$orderLG}}}{Order the SNPs in the linkage group(s).}
#' \item{\code{\link{$plotChr}}}{Plot the heatmap of the 2-point recombination fraction 
#' estimates (or LOD scores) when SNPs are ordered according to the genome assembly.}
#' \item{\code{\link{$plotLG}}}{Plot the heatmap of the 2-point recombination fraction 
#' estimates (or LOD scores) when SNPs are ordered according to their linkage groups.}
#' \item{\code{\link{$plotLM}}}{Plot linkage maps.}
#' \item{\code{\link{$plotSyn}}}{Produce a Synteny plot.}
#' \item{\code{\link{$print}}}{Print summary information for BC object.}
#' \item{\code{\link{$removeLG}}}{Remove linkage group(s).}
#' \item{\code{\link{$removeSNP}}}{Remove SNP(s) from linkage group(s).}
#' \item{\code{\link{$rf_2pt}}}{Compute the 2-point recombination fraction (and LOD score) between all SNP pairs.}
#' \item{\code{\link{$unmaskSNP}}}{Unmask SNP(s) in the data set.}
#' \item{\code{\link{$writeLM}}}{Write linkage mapping results to a file.}
#' }
#' @format NULL
#' @author Timothy P. Bilton
#' @name BC
#' @export

### R6 class for creating a data format for full-sib families
BC <- R6::R6Class("BC",
              inherit = GUSbase::RA,
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
                print = function(what = NULL, ...){
                  if(is.null(what)){
                    what <- c(what, "data")
                    if(!is.null(private$LG_mat) || !is.null(private$LG_pat))
                      what <- c(what, "LG-pts")
                    if(!is.null(private$LG_mat_bi) || !is.null(private$LG_pat_bi))
                      what <- c(what, "LG-bi")
                    if(!is.null(private$para))
                      what <- c(what, "map")
                  } else if(!is.vector(what) || !is.character(what) ||
                            any(!(what %in% c("data","LG-pts","LG-bi","map"))))
                    stop("Argument for what output to print is invalid")
                  ## print the required output
                  if(private$noFam == 1){
                    if(any(what == "data"))
                      cat(private$summaryInfo$data, sep="")
                    if(any(what == "LG-pts")){
                      if(is.null(private$LG_mat) & is.null(private$LG_pat))
                        warning("Linkage groups have not been formed. Use the '$createLG' function to form linkage groups.")
                      else{
                        cat("Pseudo-testcross Linkage Group Summary:\n")
                        MI <- unlist(lapply(private$LG_mat, length))
                        PI <- unlist(lapply(private$LG_pat, length))
                        tab <- cbind(LG=1:(length(MI)+length(PI)),MI=c(MI,rep(0,length(PI))),PI=c(rep(0,length(MI)),PI))
                        tab <- stats::addmargins(tab, margin = 1)
                        rownames(tab) <- NULL
                        tab[nrow(tab), 1] <- "TOTAL"
                        prmatrix(tab, rowlab = rep("",nrow(tab)), quote=F)
                      }
                    }
                    if(any(what == "LG-bi")){
                      if(is.null(private$LG_mat_bi) && is.null(private$LG_pat_bi))
                        warning("Linkage groups do not have any BI SNPs. Use the '$addBIsnps' to add BIsnps to the linkage groups.")
                      else{
                        cat("Pseudo-testcross Linkage Group with BI SNPs Summary:\n")
                        LGlist <- c(private$LG_mat_bi,private$LG_pat_bi)
                        ## maternal
                        if(!is.null(private$LG_mat_bi)){
                          cat("Maternal LGs:\n")
                          nmat <- 1:length(private$LG_mat_bi)
                          MI <- unlist(lapply(LGlist[nmat], function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) %in% c(4,5), na.rm=T)))
                          BI <- unlist(lapply(LGlist[nmat], function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) == 1, na.rm=T)))
                          TOTAL <- unlist(lapply(LGlist[nmat], length))
                          tab <- cbind(LG=nmat,MI,BI,TOTAL)
                          tab <- stats::addmargins(tab, margin = 1)
                          rownames(tab) <- NULL
                          tab[nrow(tab), 1] <- "TOTAL"
                          prmatrix(tab, rowlab = rep("",nrow(tab)), quote=F)
                        }
                        ## paternal
                        if(!is.null(private$LG_pat_bi)){
                          cat("Paternal LGs:\n")
                          npat <- length(private$LG_mat_bi) + 1:length(private$LG_pat_bi)
                          PI <- unlist(lapply(LGlist[npat], function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) %in% c(2,3), na.rm=T)))
                          BI <- unlist(lapply(LGlist[npat], function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) == 1, na.rm=T)))
                          TOTAL <- unlist(lapply(LGlist[npat], length))
                          tab <- cbind(LG=npat,PI,BI,TOTAL)
                          tab <- stats::addmargins(tab, margin = 1)
                          rownames(tab) <- NULL
                          tab[nrow(tab), 1] <- "TOTAL"
                          prmatrix(tab, rowlab = rep("",nrow(tab)), quote=F)
                        }
                      }
                    }
                    if(any(what == "map")){
                      if(is.null(private$para))
                        warning("no maps have been estimated. Use the '$computeMap' function to compute linkage maps.")
                      else{
                        cat(private$summaryInfo$map[[1]])
                        cat(private$summaryInfo$map[[2]])
                        prmatrix(private$summaryInfo$map[[3]], rowlab = rep("",nrow(private$summaryInfo$map[[3]])), quote=F)
                        if(length(private$summaryInfo$map) == 5){
                          cat(private$summaryInfo$map[[4]])
                          prmatrix(private$summaryInfo$map[[5]], rowlab = rep("",nrow(private$summaryInfo$map[[5]])), quote=F)
                        }
                      }
                    }
                  }
                  else
                    cat("Not yet implemented")
                  return(invisible(NULL))
                },
                #############################################################
                ## Function for removing SNPs from the linkage groups
                removeSNP = function(snps, where = NULL){
                  ## some checks
                  if( !is.vector(snps) || !is.numeric(snps)  || any(is.na(snps)) || !all(snps == round(snps)) || any(snps < 1) || any(snps > private$nSnps) )
                    stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                  snps <- unique(snps) ## make sure no double ups in SNPs
                  if(is.null(where)){
                    if(is.null(private$LG_mat_bi) & is.null(private$LG_pat_bi)) where = "LG-pts"
                    else where = "LG-comb"
                  }
                  if(where == "LG-pts"){
                    if (is.null(private$LG_mat) & is.null(private$LG_pat))
                      stop("All linkage groups have been removed. Use the '$createLG' function to create linkage groups.")
                    else{
                      ## remove SNP from the maternal linkage groups
                      for(lg in 1:length(private$LG_mat)){
                        if(any(private$LG_mat[[lg]] %in% snps)){
                          lgremove <- which(private$LG_mat[[lg]] %in% snps)
                          if(length(lgremove) > 0)
                            private$LG_mat[[lg]] <- private$LG_mat[[lg]][-lgremove]
                        }
                      }
                      ## remove SNP from the paternal linkage groups
                      for(lg in 1:length(private$LG_pat)){
                        if(any(private$LG_pat[[lg]] %in% snps)){
                          lgremove <- which(private$LG_pat[[lg]] %in% snps)
                          if(length(lgremove))
                            private$LG_pat[[lg]] <- private$LG_pat[[lg]][-lgremove]
                        }
                      }
                      ## check that there are no LGs that are now empty
                      empty_mat <- lapply(private$LG_mat, length) == 0
                      if(any(empty_mat))
                        private$LG_mat[which(empty_mat)] <- NULL
                      empty_pat <- lapply(private$LG_pat, length) == 0
                      if(any(empty_pat))
                        private$LG_pat[which(empty_pat)] <- NULL
                      ## check whether there are no more PI or MI linkage groups left
                      if(length(private$LG_mat) == 0)
                        private$LG_mat <- NULL
                      if(length(private$LG_pat) == 0)
                        private$LG_pat <- NULL
                    }
                  }
                  else if(where == "LG-bi"){
                    if (is.null(private$LG_mat_bi) & is.null(private$LG_pat_bi))
                      stop("No linakge groups with BI SNPs. Use the '$addBIsnps' function to add BI SNPs to linkage groups.")
                    else{
                      ## remove SNP from the maternal linkage groups
                      for(lg in 1:length(private$LG_mat_bi)){
                        if(any(private$LG_mat_bi[[lg]] %in% snps)){
                          lgremove <- which(private$LG_mat_bi[[lg]] %in% snps)
                          if(length(lgremove) > 0)
                            private$LG_mat_bi[[lg]] <- private$LG_mat_bi[[lg]][-lgremove]
                        }
                      }
                      ## remove SNP from the paternal linkage groups
                      for(lg in 1:length(private$LG_pat_bi)){
                        if(any(private$LG_pat_bi[[lg]] %in% snps)){
                          lgremove <- which(private$LG_pat_bi[[lg]] %in% snps)
                          if(length(lgremove))
                            private$LG_pat_bi[[lg]] <- private$LG_pat_bi[[lg]][-lgremove]
                        }
                      }
                      ## check that there are no LGs that are now empty
                      empty_mat <- lapply(private$LG_mat_bi, length) == 0
                      if(any(empty_mat))
                        private$LG_mat_bi[which(empty_mat)] <- NULL
                      empty_pat <- lapply(private$LG_pat_bi, length) == 0
                      if(any(empty_pat))
                        private$LG_pat_bi[which(empty_pat)] <- NULL
                      ## check whether there are no more PI or MI linkage groups left
                      if(length(private$LG_mat_bi) == 0)
                        private$LG_mat_bi <- NULL
                      if(length(private$LG_pat_bi) == 0)
                        private$LG_pat_bi <- NULL
                    }
                  }
                  else
                    stop("invalid second argument ('where').")
                  return(invisible(NULL))
                },
                ## Function for removing linkage groups 
                removeLG = function(LG, where = NULL){
                  ## Check the where input
                  if(is.null(where)){
                    if(is.null(private$LG_mat_bi) & is.null(private$LG_pat_bi)) where = "LG-pts"
                    else where = "LG-bi"
                  }
                  if(where == "LG-bi"){
                    if(is.null(private$LG_mat_bi) & is.null(private$LG_pat_bi))
                      stop("There are no linkage groups with BI SNPs. Use the '$addBIsnps' to create add BI SNPs to linkage groups.")
                    nmat <- length(private$LG_mat_bi)
                    npat <- length(private$LG_pat_bi)
                    if(isValue(LG, type="pos_integer", minv=1,maxv=nmat+npat))
                      stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nmat+npat))
                    else{
                      private$LG_mat_bi <- private$LG_mat[which(!(1:nmat %in% LG))]
                      private$LG_pat_bi <- private$LG_pat[which(!((nmat+1:npat) %in% LG))]
                    }
                    if(length(private$LG_mat_bi) == 0)
                      private$LG_mat <- NULL
                    if(length(private$LG_pat_bi) == 0)
                      private$LG_pat <- NULL
                  }  else if(where == "LG-pts"){
                    if(is.null(private$LG_mat) & is.null(private$LG_pat))
                      stop("No linkage groups exist. Please use the $createLG function to create some linkage groups")
                    nmat <- length(private$LG_mat)
                    npat <- length(private$LG_pat)
                    if(isValue(LG, type="pos_integer", minv=1, maxv=nmat+npat))
                      stop(paste0("The LG indices must an integer vector between 1 and the number of linkage groups which is ",nmat+npat))
                    else{
                      private$LG_mat <- private$LG_mat[which(!(1:nmat %in% LG))]
                      private$LG_pat <- private$LG_pat[which(!((nmat+1:npat) %in% LG))]
                    }
                    if(length(private$LG_mat) == 0)
                      private$LG_mat <- NULL
                    if(length(private$LG_pat) == 0)
                      private$LG_pat <- NULL
                  }
                  else
                    stop("invalid second argument ('where').")
                  return(invisible(NULL))
                },
                ## function for merging linkage groups
                mergeLG = function(LG, where = NULL, mergeTo=NULL){
                  ## checks
                  if(!is.null(mergeTo) && (!is.character(mergeTo) || length(mergeTo) != 1 || !(mergeTo %in% c("maternal","paternal","none"))))
                    stop("Argument for which parent group is to be merged to is invalid")
                  ## Check the where input
                  if(is.null(where)){
                    if(is.null(private$LG_mat_bi) && is.null(private$LG_pat_bi)) where = "LG-pts"
                    else where = "LG-bi"
                  }
                  if(where == "LG-pts"){
                    ## Check that we have linkage groups
                    if(is.null(private$LG_mat) & is.null(private$LG_pat))
                      stop("No linkage groups exist. Please use the `$createLG` function to create some linkage groups")
                    nmat <- length(private$LG_mat)
                    npat <- length(private$LG_pat)
                    ## check input
                    if(isValue(LG, type="pos_integer", minv=1, maxv=nmat+npat))
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
                        ## Added merged linkage groups
                        private$LG_mat[[matgroups[1]]] <- mergedLG
                        private$LG_mat[matgroups[-1]] <- NULL
                        ## remove paternal linkage groups and change config
                        added_pat <- private$LG_pat[patgroups][[1]][private$LG_pat[patgroups][[1]] %in% c(private$group$MI,private$group$PI)]
                        if(length(added_pat) > 0)
                          private$config[[1]][added_pat] <- c(private$config[[1]][added_pat] %% 2) + 4
                        added_pat_infer <- private$LG_pat[patgroups][[1]][private$LG_pat[patgroups][[1]] %in% private$group_infer$SI]
                        if(length(added_pat_infer) > 0)
                          private$config_infer[[1]][added_pat_infer] <- c(private$config_infer[[1]][added_pat_infer] %% 2) + 4
                        #private$group$PI <- setdiff(private$group$PI, private$LG_pat[patgroups[1]][[1]])
                        #private$group$MI <- sort(c(private$group$MI,private$LG_pat[patgroups[1]][[1]]))
                        private$LG_pat[patgroups] <- NULL
                      } else if(mergeTo == "paternal"){
                        ## add merged linkage groups
                        private$LG_pat[[patgroups[1]]] <- mergedLG
                        private$LG_pat[patgroups[-1]] <- NULL
                        ## remove maternal linkage groups and change config
                        added_mat <- private$LG_mat[matgroups][[1]][private$LG_mat[matgroups][[1]] %in% c(private$group$MI,private$group$PI)]
                        if(length(added_mat) > 0)
                          private$config[[1]][added_mat] <- c(private$config[[1]][added_mat] %% 2) + 2
                        added_mat_infer <- private$LG_mat[matgroups][[1]][private$LG_mat[matgroups][[1]] %in% private$group_infer$SI]
                        if(length(added_mat_infer) > 0)
                          private$config_infer[[1]][added_mat_infer] <- c(private$config_infer[[1]][added_mat_infer] %% 2) + 2
                        #private$group$MI <- setdiff(private$group$MI, private$LG_mat[matgroups[1]][[1]])
                        #private$group$PI <- sort(c(private$group$PI,private$LG_mat[matgroups[1]][[1]]))
                        private$LG_mat[matgroups] <- NULL
                      }
                    } else{
                      ## remove the LGs that were merged
                      if(length(matgroups) > 1)
                        private$LG_mat[matgroups[-1]] <- NULL
                      if(length(patgroups) > 1)
                        private$LG_pat[patgroups[-1]] <- NULL
                    }
                  } else if(where == "LG-bi"){
                    if(is.null(private$LG_mat_bi) & is.null(private$LG_pat_bi))
                      stop("There are no linkage groups with BI SNPs. Use the '$addBIsnps' to create linkage groups with BI snps.")
                    nmat <- length(private$LG_mat_bi)
                    npat <- length(private$LG_pat_bi)
                    ## check input
                    if(isValue(LG, type="pos_integer", minv=1,maxv=nmat+npat))
                      stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nmat+npat))
                    LG <- unique(LG)
                    if(length(LG) == 1)
                      stop("Only one unique linkage group supplied. Need at least two linkage groups to merge")
                    matgroups <- LG[which(LG %in% 1:nmat)]
                    patgroups <- LG[which(LG %in% (nmat + 1:npat))] - nmat
                    if(length(matgroups) != 0 & length(patgroups) != 0)
                      stop("Cannot merge maternal LGs with paternal LGs")
                    ## merge maternals 
                    if(length(matgroups) > 1)
                      private$LG_mat_bi[[matgroups[1]]] <- unlist(private$LG_mat_bi[matgroups])
                    ## merge paternals
                    if(length(patgroups) > 1)
                      private$LG_pat_bi[[patgroups[1]]] <- unlist(private$LG_pat_bi[patgroups])
                    ## remove the LGs that were merged
                    if(length(matgroups) > 1)
                      private$LG_mat_bi[matgroups[-1]] <- NULL
                    if(length(patgroups) > 1)
                      private$LG_pat_bi[patgroups[-1]] <- NULL
                  } 
                  else
                    stop("invalid second argument ('where').")
                  return(invisible(NULL))
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
                  return(invisible(self))
                },
                #####################################################################
                ## Function for computing the 2-point rf estimates
                rf_2pt = function(nClust = 2, err = FALSE){
                  ## do some checks
                  if(!is.numeric(nClust) || !is.vector(nClust) || length(nClust)!=1 || nClust < 0 || is.infinite(nClust) )
                    stop("Number of clusters for parallelization needs to be a positive finite integer number")
                  if(!is.vector(err) || !is.logical(err) || length(err) != 1)
                    stop("Argument specifying whether error parameters are to be estimated is invlaid.")
                  ## If there is only one family
                  if(private$noFam == 1)
                    mat <- rf_2pt_single(private$ref[[1]], private$alt[[1]],
                                         private$config_orig[[1]], private$config_infer_orig[[1]],
                                         private$group, private$group_infer,
                                         nClust, private$nInd, err=err)
                  else
                    #stop("Multiple families have yet to be implemented")
                    mat <- rf_2pt_multi(private$ref, private$alt,
                                        private$config_orig,private$group, nClust, private$noFam)
                  ## Save the results to the object
                  private$rf <- mat$rf
                  private$LOD <- mat$LOD
                  return(invisible(NULL))
                },
                ## Function for creating linkage groups
                createLG = function(parent="both", LODthres=10, nComp=10, reset=FALSE){
                  ## Do some checks
                  if(is.null(private$rf) || is.null(private$LOD))
                    stop("Recombination fractions and LOD scores have not been computed.\n  Use $rf_2pt() to compute the recombination fractions and LOD scores.")
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop(cat("parent argument is not a string of length one or is invalid:",
                         "   Please select one of the following:",
                         "     maternal: Only MI SNPs",
                         "     paternal: Only PI SNPs",
                         "     both:     MI and PI SNPs"))
                  if(!is.numeric(LODthres) || !is.vector(LODthres) || length(LODthres) != 1 || is.na(nComp) ||  LODthres < 0 || !is.finite(LODthres))
                    stop("The LOD threshold (argument 2) needs to be a finite numeric number.")
                  if(!is.numeric(nComp) || !is.vector(nComp) || length(nComp) != 1 || is.na(nComp) || nComp < 0 || !is.finite(nComp) ||
                     round(nComp) != nComp)
                    stop("The number of comparsion points (argument 3) needs to be a finite integer number.")
                  if(!is.vector(reset) || !is.logical(reset) || length(reset) != 1 || is.na(reset))
                    stop("Argument `reset` is invalid. Should be a logical value.")
                  
                  ## Reset the groups 
                  if(reset){
                    private$config <- private$config_orig
                    private$config_infer <- private$config_infer_orig
                    private$LG_mat = private$LG_pat = private$LG_mat_bi = private$LG_pat_bi = NULL
                    private$summaryInfo <- private$summaryInfo$data
                  } else if(parent == "maternal") {
                    private$LG_mat = private$LG_mat_bi = NULL
                  } else if(parent == "paternal"){
                    private$LG_pat = private$LG_pat_bi = NULL
                  } else if(parent == "both") {
                    private$LG_mat = private$LG_pat = private$LG_mat_bi = private$LG_pat_bi = NULL
                  }
                  
                  ## Create the groups
                  if(parent == "maternal" || parent == "both")
                    private$LG_mat <- createLG(private$group, private$LOD, "maternal", LODthres, nComp, masked=private$masked)
                  if(parent == "paternal" || parent == "both")
                    private$LG_pat <- createLG(private$group, private$LOD, "paternal", LODthres, nComp, masked=private$masked)
                  private$LG_mat_bi <- private$LG_pat_bi <- NULL
                  private$LG_map <- summaryInfo <- NULL
                  ## display the results
                  self$print(what = "LG-pts")
                  return(invisible(NULL))
                },
                ## function for adding the unmapped (or inferred SNPs) to the linkage groups
                addSNPs = function(LODthres=10, nComp=10){
                  
                  ## Do some checks
                  if(!is.null(private$LG_mat) & !is.null(private$LG_pat)) parent = "both"
                  else if(!is.null(private$LG_mat)) parent = "maternal"
                  else if(!is.null(private$LG_pat)) parent = "paternal"
                  else stop("No linkage groups exist. Use the `$createLG` function to create linkage groups")
                  if(!is.numeric(LODthres) || !is.vector(LODthres) || length(LODthres) != 1 || LODthres < 0 || !is.finite(LODthres))
                    stop("The LOD threshold (argument 1) needs to be a finite positive numeric number.")
                  if(!is.numeric(nComp) || !is.vector(nComp) || length(nComp) != 1 || nComp < 0 || !is.finite(nComp) ||
                     round(nComp) != nComp)
                    stop("The number of comparsion points (argument 2) needs to be a finite positive integer number.")

                  ## new LGs for adding BI snps to
                  newLGlist <- c(private$LG_mat, private$LG_pat)
                  
                  ## Find the unmapped loci
                  
                  unmapped <- sort(unlist(c(private$group$MI[which(!(private$group$MI %in% unlist(newLGlist)))],
                                            private$group$PI[which(!(private$group$PI %in% unlist(newLGlist)))],
                                            private$group_infer$SI[which(!(private$group_infer$SI %in% unlist(newLGlist)))])))
                  # check for masked SNPs
                  unmapped <- unmapped[which(!private$masked[unmapped])]
                  
                  if(length(unmapped) == 0)
                    stop("There are no SNPs remaining that are unmapped")
                  
                  ## count the number of LGs
                  nLG <- length(newLGlist)

                  ## Run algorithm for generating the linkage groups
                  noneMapped = FALSE
                  count = 0
                  added <- numeric(0)
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
                          added <- c(added, snp)
                          newLG <- which.max(LODvalue)
                          newLGlist[[newLG]] <- c(newLGlist[[newLG]], snp)
                          unmapped <- unmapped[-which(unmapped == snp)]
                          noneMapped = FALSE
                        }
                      }
                    }
                  }
                  ## Check to see if any SNPs were added
                  if(count == 0){
                    stop("No SNPs were added to the Linkage Groups")
                  } else{
                    if(length(private$LG_mat) > 0){
                      added_mat <- added[which((added %in% unlist(c(private$group$MI,private$group$PI))) & 
                                                 (added %in% unlist(newLGlist[1:length(private$LG_mat)])))]
                    } else added_mat <- numeric(0)
                    if(length(private$LG_pat) > 0){
                      added_pat <- added[which((added %in% unlist(c(private$group$MI,private$group$PI))) & 
                                                 (added %in% unlist(newLGlist[length(private$LG_mat) + 1:length(private$LG_pat)])))]
                    } else added_pat <- numeric(0)
                    if(length(private$LG_mat) > 0){
                      added_mat_infer <- added[which( (added %in% unlist(newLGlist[1:length(private$LG_mat)])) &
                                                        (added %in% private$group_infer$SI))]
                    } else added_mat_infer <- numeric(0)
                    if(length(private$LG_pat) > 0){
                      added_pat_infer <- added[which( (added %in% unlist(newLGlist[length(private$LG_mat) + 1:length(private$LG_pat)])) &
                                                        (added %in% private$group_infer$SI))]
                    } else added_pat_infer <- numeric(0)
                    ## update configations if required
                    if(length(added_mat) > 0)
                      private$config[[1]][added_mat] <- (c(private$config[[1]][added_mat]) %% 2) + 4
                    if(length(added_pat) > 0)
                      private$config[[1]][added_pat] <- (c(private$config[[1]][added_pat]) %% 2) + 2
                    if(length(added_mat_infer) > 0)
                      private$config_infer[[1]][added_mat_infer] <- (c(private$config_infer[[1]][added_mat_infer]) %% 2) + 4
                    if(length(added_pat_infer) > 0)
                      private$config_infer[[1]][added_pat_infer] <- (c(private$config_infer[[1]][added_pat_infer]) %% 2) + 2
                    ## set the new LGs
                    if(!is.null(private$LG_mat))
                      private$LG_mat <- newLGlist[1:length(private$LG_mat)]
                    if(!is.null(private$LG_pat))
                      private$LG_pat <- newLGlist[length(private$LG_mat) + 1:length(private$LG_pat)]
                  }
                  self$print(what = "LG-pts")
                  return(invisible())
                },
                ## Function for adding the informative SNPs to the LGs
                addBIsnps = function(LODthres=10, nComp=10){
                  
                  ## Do some checks
                  if(!is.null(private$LG_mat) & !is.null(private$LG_pat)) parent = "both"
                  else if(!is.null(private$LG_mat)) parent = "maternal"
                  else if(!is.null(private$LG_pat)) parent = "paternal"
                  else stop("No linkage groups exist. Use the `$createLG` function to create linkage groups")
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
                  if(parent == "maternal" || parent == "both" ){
                    private$LG_mat_bi <- private$mapBISnps(unmapped, parent="maternal", LODthres=LODthres, nComp=nComp)
                  }
                  
                  ## if paternal groups, map the BIs
                  if(parent == "paternal" || parent == "both" ){
                    private$LG_pat_bi <- private$mapBISnps(unmapped, parent="paternal", LODthres=LODthres, nComp=nComp)
                  }
                  
                  private$para <- NULL            ## reset the computed maps
                  private$summaryInfo$map <- NULL ## reset the computed maps summary
                  self$print(what = "LG-bi")
                  return(invisible(NULL))
                },
                ## Function for ordering linkage groups
                orderLG = function(LG = NULL, mapfun = "haldane", weight="LOD2", ndim=30, spar=NULL, filename=NULL){
                  ## do some checks
                  if(is.null(private$LG_mat_bi) & is.null(private$LG_pat_bi)){
                    warning("Linkage groups do not contain BI SNPs. Proceeding with only partially informative SNPs")
                    if(!is.null(private$LG_mat) & !is.null(private$LG_pat)) parent = "both"
                    else if(!is.null(private$LG_mat)) parent = "maternal"
                    else if(!is.null(private$LG_pat)) parent = "paternal"
                    else stop("No linkage groups exist. Use the `$createLG` function to create linkage groups")
                    private$LG_mat_bi <- private$LG_mat
                    private$LG_pat_bi <- private$LG_pat
                  } else{
                    if(!is.null(private$LG_mat_bi) & !is.null(private$LG_pat_bi)) parent = "both"
                    else if(!is.null(private$LG_mat_bi)) parent = "maternal"
                    else if(!is.null(private$LG_pat_bi)) parent = "paternal"
                  }
                  LGlist <- c(private$LG_mat_bi,private$LG_pat_bi)
                  nmat <- length(private$LG_mat_bi)
                  npat <- length(private$LG_pat_bi)
                  ## Work out if we want a subset of the LGs
                  if(!is.null(LG)){
                    if(isValue(LG, type="pos_integer", minv=1, maxv=length(LGlist)))
                      stop(paste0("At least one linkage group number does not exist. Indices must be between 1 and",length(LGlist),"\n"))
                  } else LG <- 1:length(LGlist)
                  if(!(mapfun %in% c("morgan","haldane","kosambi")))
                    stop("Unknown mapping function")
                  if(!(weight %in% c("LOD","LOD2","none")))
                    stop("Unknown weighting function")
                  if(!is.null(filename) && (!is.vector(filename) || !is.character(filename) || length(filename) != 1))
                    stop("Specified filename is invalid")
                  if(!is.null(filename)){
                    temp <- file.create(paste0(filename,"_LG",LG[1],".pdf"))
                    if(!temp)
                      stop("Error in filename. Cannot create pdf")
                  }
                  if(nmat > 0) matgroups <- LG[which(LG %in% 1:nmat)]
                  else matgroups <- NULL
                  if(npat > 0) patgroups <- LG[which(LG %in% (nmat + 1:npat))]
                  else patgroups <- NULL
                  ## order the linkage groups
                  for(lg in LG){
                    ind <- LGlist[[lg]]
                    ## set up the weighting matrix
                    if(weight == "LOD")
                      wmat <- matrix(private$LOD[ind,ind],nrow=length(ind), ncol=length(ind))
                    else if (weight == "LOD2")
                      wmat <- matrix((private$LOD[ind,ind])^2,nrow=length(ind), ncol=length(ind))
                    else if (weight == "none")
                      wmat <- matrix(1,nrow=length(ind), ncol=length(ind))
                    ## set of the distance matrix
                    dmat <- mfun(private$rf[ind,ind], fun = mapfun, centiM=FALSE)
                    dmat[which(is.infinite(dmat))] <- 18.3684
                    ## unconstrainted MDS
                    MDS <- smacof::smacofSym(delta=dmat, ndim=ndim, weightmat = wmat, itmax = 1e+05)
                    ## principal curve
                    pcurve <- princurve::principal_curve(MDS$conf, maxit = 150, spar = spar)
                    ## plot the results
                    if(is.null(filename)) temp_par <- par(no.readonly=T)
                    else grDevices::pdf(paste0(filename,"_LG",lg,".pdf"), width=8, height=8)
                    graphics::par(mfrow = c(1,2))
                    graphics::plot(MDS$conf[,1],MDS$conf[,2], ylab="Dimension 2", xlab="Dimension 1", type="n")
                    graphics::text(MDS$conf, labels=ind)
                    graphics::lines(pcurve)
                    graphics::image(private$rf[ind,ind][pcurve$ord,pcurve$ord], axes=F, 
                                    col=grDevices::heat.colors(100))
                    if(is.null(filename)) graphics::par(temp_par)
                    else grDevices::dev.off()
                    ## Set the new order
                    if(lg %in% matgroups)
                      private$LG_mat_bi[[lg]] <- ind[pcurve$ord]
                    else if (lg %in% patgroups)
                      private$LG_pat_bi[[lg - nmat]] <- ind[pcurve$ord]
                  }
                  ## Acknowledge paper we are using algorithm from
                  cat("Using the MDS scaling algorithm to order SNPs.\n",
                      "Please cite:\n",
                      "\tPreedy & Hackett (2016). Theor Appl Genet. 129:2117-2132\n",sep="")
                  
                  return(invisible(NULL))
                },
                ## Function for plotting linkage groups
                plotLG = function(parent = "both", LG=NULL, mat="rf", filename=NULL, interactive=FALSE, what = NULL, ...){
                  ## do some checks
                  if(!is.vector(mat) || !is.character(mat) || length(mat) != 1 || !(mat %in% c('rf','LOD')))
                    stop("Argument specifying which matrix to plot (argument 1) must be either 'rf' or 'LOD'")
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop(paste("Parent argument is not a string of length one or is invalid:",
                               "   Please select one of the following:",
                               "   maternal: Only MI and BI SNPs",
                               "   paternal: Only PI and BI SNPs",
                               "   both:     MI, PI and BI SNPs", sep="\n"))
                  if(!is.vector(interactive) || length(interactive) != 1 || !is.logical(interactive))
                    stop("Argument sepcifting whether to produce an interactive heatmap or not is invalid.")
                  if(!is.null(filename) && (!is.vector(filename) || !is.character(filename) || length(filename) != 1))
                    stop("Specified filename is invalid")
                  ## Check the what input
                  if(is.null(what)){
                    if(is.null(private$LG_mat_bi) && is.null(private$LG_pat_bi)) what = "LG-pts"
                    else what = "LG-bi"
                  }
                  plotly.arg <- list(...)
                  if(is.null(plotly.arg$selfcontained))
                    plotly.arg$selfcontained = FALSE
                  
                  if(private$noFam == 1){
                    ## Work out which LGs list to use
                    if(what == "LG-pts"){
                      if(is.null(private$LG_pat) && is.null(private$LG_mat))
                        stop("No linkage groups are available to be plotted. Please use the `$createLG` function to create some linkage groups")
                      else{
                        if(parent == "maternal"){
                          LGlist = private$LG_mat
                          if(length(LGlist) == 0)
                            stop("No maternal linkage groups to plot")
                          LGnum = 1:length(private$LG_mat)
                        } else if(parent == "paternal"){
                          LGlist = private$LG_pat
                          if(length(LGlist) == 0)
                            stop("No paternal linkage groups to plot")
                          LGnum = length(private$LG_mat) + 1:length(private$LG_pat)
                        } else if (parent == "both"){
                          LGlist <- c(private$LG_mat,private$LG_pat)
                          LGnum = 1:length(LGlist)
                        }
                      }
                    } else if(what == "LG-bi"){
                      if(is.null(private$LG_mat_bi) && is.null(private$LG_pat_bi))
                        stop("There are no combined linkage groups with BI SNPs. Use the '$addBIsnps' to create combined linkage groups.")
                      else{
                        if(parent == "maternal"){
                          LGlist = private$LG_mat_bi
                          if(length(LGlist) == 0)
                            stop("No maternal linkage groups to plot")
                          LGnum = 1:length(private$LG_mat_bi)
                        } else if(parent == "paternal"){
                          LGlist = private$LG_pat_bi
                          if(length(LGlist) == 0)
                            stop("No paternal linkage groups to plot")
                          LGnum = length(private$LG_mat_bi) + 1:length(private$LG_pat_bi)
                        } else if (parent == "both"){
                          LGlist <- c(private$LG_mat_bi,private$LG_pat_bi)
                          LGnum = 1:length(LGlist)
                        }
                      }
                    } else
                      stop("invalid argument 'what'")
                    ## Work out if we want a subset of the LGs
                    if(!is.null(LG)){
                      if(GUSbase::checkVector(LG, type="pos_integer", minv=min(LGnum), maxv=max(LGnum)))
                        stop(paste0("At least one linkage group number does not exist. Indices must be between ",
                                    min(LGnum)," and ",max(LGnum)," for ",parent," LGs\n"))
                      LGlist <- LGlist[which(LGnum %in% LG)]
                      names(LGlist) <- LG
                    }
                    else
                      names(LGlist) <- LGnum
                    
                    ## Sort out the matrix
                    if(mat == "rf"){
                      temprf <- private$rf
                      mycol <- heat.colors(100)
                      zlim = c(0,0.5)
                    }
                    else if(mat == "LOD"){
                      temprf <- private$LOD
                      mycol <- grDevices::colorRampPalette(rev(c("red","orange","yellow","white")), bias=5)(100)
                      zlim=c(0,50)
                    }
                    b <- ncol(temprf) + 1
                    temprf <- cbind(temprf,rep(NA,b-1))
                    temprf <- rbind(temprf,rep(NA,b))
                    chrom.ind <- unlist(lapply(LGlist, function(x) c(x,b)), use.names = FALSE)
                    chrom.ind <- chrom.ind[-length(chrom.ind)]
                    temprf <- temprf[chrom.ind, chrom.ind]
                    nn <- length(chrom.ind)
                    b_indx <- chrom.ind == b
                    
                    if(interactive){
                      ## Plot the matrix
                      chrom.ind[which(!b_indx)] <- paste0(chrom.ind[which(!b_indx)]," (", rep(names(LGlist), lapply(LGlist,length)),")") 
                      chrom.ind[which(b_indx)] <- rep("Break",length(LGlist)-1)
                      if(mat == "rf")
                        hovertext <- matrix(paste(matrix(paste0("row: ",chrom.ind), nrow=nn, ncol=nn), 
                                                  matrix(paste0("col: ",chrom.ind), nrow=nn, ncol=nn, byrow=T), paste0("rf: ",round(temprf,4)), sep="<br>"),
                                            nrow=nn, ncol=nn)
                      else if(mat == "LOD")
                        hovertext <- matrix(paste(matrix(paste0("row: ",chrom.ind), nrow=nn, ncol=nn), 
                                                  matrix(paste0("col: ",chrom.ind), nrow=nn, ncol=nn, byrow=T), paste0("LOD: ",round(temprf,4)), sep="<br>"),
                                            nrow=nn, ncol=nn)
                      
                        
                      ax <- list(visible=FALSE)
                      ax_z <- list(range=zlim)
                      # suppress warnings  
                      storeWarn <- getOption("warn")
                      options(warn = -1) 
                      ## produce the plotly plot
                      if(length(which(b_indx)) == 0){
                        p <- plotly::plot_ly(z=temprf, type="heatmap", showscale=F, hoverinfo="text",
                                             text=hovertext, colors=mycol) %>%
                          plotly::layout(margin=list(l=0,r=0,t=0,b=0), xaxis=ax, yaxis=ax, zaxis=ax_z)
                      }
                      else{
                        p <- plotly::plot_ly(z=temprf, type="heatmap", showscale=F, hoverinfo="text",
                                             text=hovertext, colors=mycol) %>% 
                          plotly::add_segments(x=which(b_indx)-1,xend=which(b_indx)-1,y=0,yend=nn, line=list(color="black"),  showlegend=F) %>%
                          plotly::add_segments(y=which(b_indx)-1,yend=which(b_indx)-1,x=0,xend=nn, line=list(color="black"),  showlegend=F) %>%
                          plotly::layout(margin=list(l=0,r=0,t=0,b=0), xaxis=ax, yaxis=ax, zaxis=ax_z)
                      }
                      if(!is.null(filename)){
                        temp <- file.create(paste0(filename,".html"))
                        if(!temp)
                          stop("Error in filename. Cannot create png")
                        else{
                          htmlwidgets::saveWidget(p, paste0(filename,".html"), selfcontained = plotly.arg$selfcontained)
                        }  
                      }
                      else 
                        print(p)
                      options(warn = storeWarn) 
                    }
                    else{
                      
                      if(!is.null(filename)){
                        temp <- file.create(paste0(filename,".png"))
                        if(!temp)
                          stop("Error in filename. Cannot create png")
                        else
                          grDevices::png(paste0(filename,".png"), width=nn,height=nn,res=72)
                      }
                      else
                        temp_par <- graphics::par(no.readonly = TRUE)
                      graphics::par(mar=rep(0,4))
                      #graphics::par(mar=rep(0,4), xaxt='n',yaxt='n',bty='n',ann=F)
                      graphics::image(x=1:nn, y=1:nn, z=temprf, zlim=zlim, col=mycol, axes=F, xlim=c(0,nn), ylim=c(1,nn+1))
                      graphics::abline(v=which(b_indx))
                      graphics::abline(h=which(b_indx))
                      if(!is.null(filename))
                        dev.off()
                      else
                        graphics::par(temp_par)
                    }
                  }
                  else
                    stop("Not yet implemented.")
                  return(invisible())
                },
                ## Function for plotting chromosome in their original ordering
                plotChr = function(chrom=NULL, parent = "maternal", mat=c("rf"), filename=NULL, chrS=2, lmai=0.25){
                  ## do some checks
                  if(!is.vector(mat) || !is.character(mat) || length(mat) != 1 || !(mat %in% c('rf','LOD')))
                    stop("Argument specifying which matrix to plot (argument 1) must be either 'rf' or 'LOD'")
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop(paste("Parent argument is not a string of length one or is invalid:",
                               "   Please select one of the following:",
                               "   maternal: Only MI and BI SNPs",
                               "   paternal: Only PI and BI SNPs",
                               "   both:     MI, PI and BI SNPs", sep="\n"))
                  nChrom = length(unique(private$chrom))
                  if(is.null(chrom)){
                    chrom <- 1:nChrom
                  } else if(GUSbase::checkVector(chrom, type = "pos_integer", maxv=nChrom))
                    stop(paste0("Invalid chromosome number (first argument). Must be an integer number between 0 and ",nChrom,".\n"))
                  
                  if(is.null(filename))
                    temp_par <- par(no.readonly = TRUE) # save the current plot margins
                  
                  ## workout which SNPs to plot to plot
                  if(private$noFam == 1){
                    ## workout the indices for the chromosomes
                    names <- unique(private$chrom)
                    if(parent == "maternal")
                      LGlist <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,4,5))), simplify=F)
                    else if(parent == "paternal")
                      LGlist <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,2,3))), simplify=F)
                    else if(parent == "both")
                      LGlist <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,2,3,4,5))), simplify=F)
                    ## plot the chromsomes rf info
                    if(mat == "rf")
                      plotLG(mat=private$rf, LG=LGlist[chrom], filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=T, type="rf")
                    else if(mat == "LOD")
                      plotLG(mat=private$LOD, LG=LGlist[chrom], filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=T, type="LOD")
                    else
                      stop("Matrix to be plotted not found.") ## shouldn't get here
                    if(is.null(filename))
                      graphics::par(temp_par) # reset the plot margins
                    return(invisible())
                  }
                  else{
                    stop("not implemented yet")
                  }
                },
                ## Function for plotting the linkage groups
                plotLM = function(LG = NULL, fun="haldane", col="black"){
                  if(is.null(private$LG_map))
                    stop("There are no linkage maps availabe to plot. Use the '$computeMap' function to compute linkage maps.")
                  nLGs = length(private$LG_map)
                  if(is.null(LG))
                    LG <- 1:nLGs
                  if(GUSbase::checkVector(LG, type="pos_integer", minv=1, maxv=nLGs))
                    stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nLGs))
                  matlg <- which(unlist(lapply(private$para$OPGP[LG], function(x) all(x %in% c(1:4,9:12)))))
                  patlg <- which(unlist(lapply(private$para$OPGP[LG], function(x) all(x %in% c(1:8)))))
                  LG = c(matlg,patlg)
                  nLGs = length(LG)
                  if(nLGs == 0)
                    stop("There are no linkage maps for the specified linkage groups.")
                  if(!is.vector(fun) || length(fun) != 1 || any(!(fun %in% c("morgan", "haldane", "kosambi"))))
                    stop("Input argument for mapping distance is invalid")
                  ## check the color input
                  if(!is.vector(is.character(col)) || !is.character(col))
                    stop("Input argument for colors of the linkage maps is invalid.")
                  else{
                    tc <- unname(sapply(col, function(y) tryCatch(is.matrix(grDevices::col2rgb(y)), error = function(e) FALSE)))
                    if(any(!tc))
                      stop("At least one color in the color list is undefined")
                    else
                      col = grDevices::rgb(t(grDevices::col2rgb(col)), max=255)
                    if(length(col) == 1)
                      col <- rep(col,nLGs)
                    else if(length(col) != nLGs)
                      stop("Number of colors specified does not equal the number of linkage groups")
                  }
                  temp_par <- graphics::par(no.readonly = TRUE) # save the current plot margins
                  
                  ellipseEq_pos <- function(x) c2 + sqrt(b^2*(1-round((x-c1)^2/a^2,7)))
                  ellipseEq_neg <- function(x) c2 - sqrt(b^2*(1-round((x-c1)^2/a^2,7)))
                  mapDist = lapply(c(private$para$rf[LG]),mfun, fun=fun, centi=TRUE)
                  yCoor = max(unlist(lapply(mapDist,sum)))
                  graphics::par(mar=rep(2,4), mfrow=c(1,1))
                  graphics::plot(NULL,NULL, ylim=c(yCoor+0.01*yCoor,-0.01*yCoor),xlim=c(0,0.25*(nLGs+1)), xaxt='n',bty='n',xlab="",ylab="")
                  err = 0.05
                  ## plot each map for each linkage group
                  for(lg in 1:nLGs){
                    cent = 0.25*lg
                    graphics::segments(cent-err,c(0,cumsum(mapDist[[lg]])),cent+err,c(0,cumsum(mapDist[[lg]])),col=col[lg])
                    graphics::segments(cent-err,-0.01*yCoor,cent-err,sum(mapDist[[lg]])+0.01*yCoor,col=col[lg])
                    graphics::segments(cent+err,-0.01*yCoor,cent+err,sum(mapDist[[lg]])+0.01*yCoor,col=col[lg])
                    ## Plot the curves at the ends
                    a=err; b=0.02*yCoor; c1=cent; c2=-0.01*yCoor;
                    graphics::curve(ellipseEq_neg, from=cent-err,to=cent+err,add=T,col=col[lg])
                    c2=sum(mapDist[[lg]])+0.01*yCoor;
                    graphics::curve(ellipseEq_pos, from=cent-err,to=cent+err,add=T,col=col[lg])
                  }
                  graphics::par(temp_par) # reset the plot margins
                  #if(!is.null(names))
                  #  mtext(names, side=1, at=seq(0.25,0.25*nLGs,0.25))
                },
                ## compare order with original and ordered markers
                plotSyn = function(){
                  if(is.null(private$LG_mat_bi) && is.null(private$LG_pat_bi))
                    stop("There are no ordered linkage groups availabe. Use the '$orderLG' function to order linkage groups.")
                  ## work out the LG and original orders
                  LGlist <- c(private$LG_mat_bi,private$LG_pat_bi)
                  LGorder <- unlist(LGlist)
                  orgOrder <- sort(LGorder)
                  ## determine the breask on the plot
                  temp <- cumsum(unlist(lapply(LGlist, length)))
                  LGbreaks <- orgOrder[temp[-length(temp)]] + 0.5
                  chrBreaks <- which(diff(as.numeric(as.factor(private$chrom)))==1) + 0.5
                  ## plot the synteny plot
                  temp_par <- graphics::par(no.readonly = TRUE) # save the current plot margins
                  graphics::par(mfrow=c(1,1))
                  graphics::plot(orgOrder,LGorder, pch=20,cex=0.8, xaxt="n", yaxt="n",ylab="Assembly Order", xlab="Linkage Group Order", 
                                 ylim=c(min(orgOrder),max(orgOrder)), xlim=c(min(LGorder), max(LGorder)), bty='n')
                  graphics::abline(v=c(min(LGorder)-1,LGbreaks,max(LGorder)+1))
                  graphics::abline(h=c(min(orgOrder)-1,chrBreaks,max(orgOrder)+1))
                  graphics::abline(v=LGbreaks)
                  graphics::mtext(text = unique(private$chrom[orgOrder]), side = 2, 
                                  at = apply(cbind(c(min(orgOrder),chrBreaks),c(chrBreaks,max(orgOrder))),1,mean))
                  graphics::mtext(text = 1:length(LGlist), side = 1, 
                                  at = apply(cbind(c(min(LGorder),LGbreaks),c(LGbreaks,max(LGorder))),1,mean))
                  graphics::par(temp_par) # reset the plot margins
                }, 
                ## Function for computing the rf's for each chromosome 
                computeMap = function(parent = "both", chrom=NULL, init_r=0.001, ep=0.001, method=NULL, err=TRUE, multiErr=FALSE, 
                                      mapped=TRUE, nThreads=1, inferOPGP=TRUE, rfthres=0.1){
                  ## do some checks
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop(paste("parent argument is not a string of length one or is invalid:",
                              "   Please select one of the following:",
                              "   maternal: Only MI and BI SNPs",
                              "   paternal: Only PI and BI SNPs",
                              "   both:     MI, PI and BI SNPs", sep="\n"))
                  if( !is.null(init_r) & (GUSbase::checkVector(init_r, type="pos_numeric", minv = 0, maxv=0.4) || length(init_r) != 1))
                    stop("Starting values for the recombination fraction needs to be a numeric vector or integer or a NULL object")
                  if( ((length(ep) != 1)  || !is.numeric(ep) || (ep <= 0 | ep >= 1)) )
                    stop("Value for the error parameters needs to be a single numeric value in the interval (0,1) or a NULL object")
                  if(!is.null(method) && (!is.vector(method) || !is.character(method) || 
                                          length(method) != 1 || !(method %in% c("optim", "EM"))))
                    stop("Argument `method` is invalid. Must be NULL, 'optim' or 'EM'")
                  if(!is.logical(err) || length(err) != 1 || is.na(err))
                    stop("Argument `err` is invalid. Must be a logical value")
                  if(!is.logical(multiErr) || length(multiErr) != 1 || is.na(multiErr))
                    stop("Argument `multiErr` is invalid. Must be a logical value")
                  if(!is.logical(inferOPGP) || length(inferOPGP) != 1 || is.na(inferOPGP))
                    stop("Argument `inferOPGP` is invalid. Must be a logical value")                  
                  if(!is.logical(mapped) || length(mapped) != 1 || is.na(mapped))
                    stop("Argument `mapped` is invalid. Must be a logical value")                  
                  if(GUSbase::checkVector(nThreads, type="pos_integer", minv=1) || length(nThreads) != 1)
                    stop("Argument `nThreads` is invalid. Must be a positive integer")
                  if(GUSbase::checkVector(rfthres, type="pos_numeric", minv=0, maxv = 0.5) || length(rfthres) != 1)
                    stop("Argument `rfthres` is invalid. Must be a positive integer")
                  
                  ## for existing chromosome orders
                  if(!mapped){
                    stop("yet to be implemented")
                    # if(!is.null(chrom) && !is.vector(chrom))
                    #   stop("chromosomes names must be a character vector.")
                    # if(is.null(chrom)) # compute distances for all chromosomes
                    #   chrom <- unique(private$chrom)
                    # else if(!(any(chrom %in% private$chrom[!private$masked])))
                    #   stop("At least one chromosome not found in the data set.")
                    # else
                    #   chrom <- as.character(chrom) # ensure that the chromsome names are characters
                    # nchrom <- length(unique(private$chrom))
                    # cat("Computing recombination fractions:\n")
                    # if(is.null(private$para))
                    #   private$para <- list(OPGP=vector(mode = "list",length = nchrom),rf_p=vector(mode = "list",length = nchrom),
                    #                        rf_m=vector(mode = "list",length = nchrom),ep=vector(mode = "list",length = nchrom),
                    #                        loglik=vector(mode = "list",length = nchrom))
                    # else if(length(private$para$rf_p) != nchrom)
                    #   private$para <- list(OPGP=vector(mode = "list",length = nchrom),rf_p=vector(mode = "list",length = nchrom),
                    #                        rf_m=vector(mode = "list",length = nchrom),ep=vector(mode = "list",length = nchrom),
                    #                        loglik=vector(mode = "list",length = nchrom))
                    # for(i in chrom){
                    #   cat("chromosome: ",i,"\n")
                    #   indx_chrom <- which((private$chrom == i) & !private$masked)
                    #   ref_temp <- lapply(private$ref, function(x) x[,indx_chrom])
                    #   alt_temp <- lapply(private$alt, function(x) x[,indx_chrom])
                    #   ## estimate OPGP's
                    #   curOPGP <- private$para$OPGP[i]
                    #   if(inferOPGP || !is.list(curOPGP) || length(curOPGP[[1]]) != length(indx_chrom)){
                    #     tempOPGP <- list()
                    #     for(fam in 1:private$noFam){
                    #       tempOPGP <- c(tempOPGP,list(as.integer(infer_OPGP_FS(ref_temp[[fam]],alt_temp[[fam]],private$config[[fam]][indx_chrom], method=method, nThreads=nThreads))))
                    #     }
                    #     private$para$OPGP[i] <- tempOPGP
                    #   }
                    #   ## estimate the rf's
                    #   MLE <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[1],
                    #                    sexSpec=sexSpec, seqErr=err, method=method, nThreads=nThreads, multiErr=multiErr)
                    #   if(sexSpec){
                    #     private$para$rf_p[i]   <- list(MLE$rf_p)
                    #     private$para$rf_m[i]   <- list(MLE$rf_m)
                    #   } else{
                    #     private$para$rf_p[i]   <- list(MLE$rf)
                    #     private$para$rf_m[i]   <- list(MLE$rf)
                    #   }
                    #   private$para$ep[i]    <- list(list(MLE$ep))
                    #   private$para$loglik[i]<- list(MLE$loglik)
                    # }
                  }
                  else{
                    LGlist <- c(private$LG_mat_bi,private$LG_pat_bi)
                    if(!is.null(private$LG_mat_bi)) matlg <- 1:length(private$LG_mat_bi)
                    else matlg <- NULL
                    if(!is.null(private$LG_pat_bi)) patlg <- length(private$LG_mat_bi) + 1:length(private$LG_pat_bi)
                    else patlg <- NULL
                    nmatlg = length(matlg)
                    npatlg = length(patlg)
                    ## Check the input
                    if((length(LGlist) == 0))
                      stop("Linkage groups have not been ordered. Use the $orderLG function to order SNPs")
                    if(is.null(chrom)){ # compute distances for all chromosomes
                      if(parent == "maternal") chrom = matlg
                      else if (parent == "paternal") chrom = patlg 
                      else if (parent == "both") chrom = 1:length(LGlist)
                    } else if(parent == "both"){
                      if(GUSbase::checkVector(chrom, type="pos_integer", minv=1, maxv = length(LGlist)))
                        stop(paste0("Linkage group input must be an integer vector between 1 and ",length(LGlist),
                                    " (the total number of linkage groups)"))
                    } else if (parent == "maternal"){
                      if(GUSbase::checkVector(chrom, type="pos_integer", minv=min(matlg), maxv = max(matlg)))
                        stop(paste0("Linkage group input must be an integer vector between ",min(matlg)," and ", max(matlg),
                                    " for maternal linkage groups"))
                    } else if (parent == "paternal"){
                      if(GUSbase::checkVector(chrom, type="pos_integer", minv=min(patlg), maxv = max(patlg)))
                        stop(paste0("Linkage group input must be an integer vector between ",min(patlg)," and ", max(patlg),
                                    " for paternal linkage groups"))
                    }
                    
                    ## Compute rf's
                    cat("Computing recombination fractions:\n")
                    if(is.null(private$para) || any(lapply(private$para, length) != length(LGlist)))
                      private$para <- list(OPGP=replicate(length(LGlist), NA, simplify=FALSE),
                                           rf=replicate(length(LGlist), NA, simplify=FALSE),
                                           ep=replicate(length(LGlist), NA, simplify=FALSE),
                                           loglik=replicate(length(LGlist), NA, simplify=FALSE))
                    #else if(length(private$para$rf_p) != length(LGlist))
                    if(is.null(private$LG_map) || length(private$LG_map) != length(LGlist))
                      private$LG_map <- vector(mode="list", length = length(LGlist))
                    for(i in chrom){
                      cat("Linkage Group: ",i,"\n")
                      indx_chrom <- LGlist[[i]]
                      ref_temp <- lapply(private$ref, function(x) x[,indx_chrom])
                      alt_temp <- lapply(private$alt, function(x) x[,indx_chrom])
                      config_temp <- lapply(private$config, function(x) x[indx_chrom])
                      for(fam in 1:private$noFam){
                        tempIndx <- !is.na(private$config_infer[[fam]][indx_chrom])
                        config_temp[[fam]][tempIndx] <- private$config_infer[[fam]][indx_chrom[tempIndx]]
                      }
                      
                      ## if maternal linkage group
                      ## estimate OPGP's
                      curOPGP <- private$para$OPGP[i]
                      if(inferOPGP || !is.list(curOPGP) || length(curOPGP[[1]]) != length(indx_chrom)){
                        tempOPGP <- list()
                        for(fam in 1:private$noFam){
                          tempOPGP <- c(tempOPGP,list(as.integer(infer_OPGP_FS(ref_temp[[fam]],alt_temp[[fam]],config_temp[[fam]], method="EM", nThreads=nThreads))))
                        }
                        private$para$OPGP[i] <- tempOPGP
                      }
                      ## estimate the rf's
                      if(!is.null(method))
                        MLE <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                         sexSpec=F, seqErr=err, method=method, nThreads=nThreads, multiErr=multiErr)
                      else{
                        MLE_init <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                              sexSpec=F, seqErr=err, method="EM", nThreads=nThreads, multiErr=multiErr, maxit=20)
                        MLE <- rf_est_FS(init_r=MLE_init$rf, ep=MLE_init$ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                         sexSpec=F, seqErr=err, method="optim", nThreads=nThreads, multiErr=multiErr)
                      }
                      
                      private$para$rf[i]   <- list(MLE$rf)
                      private$para$ep[i]       <- list(MLE$ep)
                      private$para$loglik[i]   <- list(MLE$loglik)
                      private$LG_map[[i]]      <- LGlist[[i]]
                      ### Check for SNPs with really large rf values
                      highrf = which(MLE$rf > rfthres)
                      if(length(highrf) > 0){
                        snp1 = indx_chrom[highrf]
                        snp2 = indx_chrom[highrf+1]
                        junk = sapply(1:length(highrf),function(x) {
                          cat("RF-", highrf[x]," is ", round(MLE$rf[highrf[x]],4)," (between SNPs ",snp1[x]," and ",snp2[x],")\n", sep="")
                          return(invisible())
                        })
                      }
                    }
                  }
                  ## Create summary of results:
                  templist <- list("Linkage Map Summaries:\n")
                  ## maternal
                  if(!is.null(matlg)){
                    MI <- unlist(lapply(private$LG_map[matlg], function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) %in% c(4,5), na.rm=T)))
                    BI <- unlist(lapply(private$LG_map[matlg], function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) == 1, na.rm=T)))
                    TOTAL <- unlist(lapply(private$LG_map[matlg], length))
                    tab1 <- cbind(LG=matlg,MI,BI,TOTAL)
                    DIST <- extendVec(unlist(lapply(private$para$rf[matlg], function(x) round(sum(mfun(x, fun="haldane", centiM = T)),2))), nrow(tab1))
                    ERR <- extendVec(round(unlist(lapply(private$para$ep[matlg],mean)),5), nrow(tab1))
                    tab1 <- cbind(tab1,DIST,ERR)
                    tab1 <- stats::addmargins(tab1, margin = 1)
                    rownames(tab1) <- NULL
                    tab1[nrow(tab1), 1] <- "TOTAL"
                    tab1[nrow(tab1),6] <- ""
                    templist <- c(templist,list("\nMaternal LGs:\n"), list(tab1))
                  }
                  ## paternal
                  if(!is.null(patlg)){
                    PI <- unlist(lapply(private$LG_map[patlg], function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) %in% c(2,3), na.rm=T)))
                    BI <- unlist(lapply(private$LG_map[patlg], function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) == 1, na.rm=T)))
                    TOTAL <- unlist(lapply(private$LG_map[patlg], length))
                    tab2 <- cbind(LG=patlg,PI,BI,TOTAL)
                    DIST <- extendVec(unlist(lapply(private$para$rf[patlg], function(x) round(sum(mfun(x, fun="haldane", centiM = T)),2))), nrow(tab2))
                    ERR <- extendVec(round(unlist(lapply(private$para$ep[patlg],mean)),5), nrow(tab2))
                    tab2 <- cbind(tab2,DIST,ERR)
                    tab2 <- stats::addmargins(tab2, margin = 1)
                    rownames(tab2) <- NULL
                    tab2[nrow(tab2), 1] <- "TOTAL"
                    tab2[nrow(tab2),6] <- ""
                    templist <- c(templist,list("\nPaternal LGs:\n"), list(tab2))
                  }
                  private$summaryInfo$map <- templist
                  self$print(what = "map")
                  return(invisible(NULL))
                },
                ## Ratio of alleles for heterozygous genotype calls (observed vs expected)
                # Slightly different than the one in RA class
                cometPlot = function(filename=NULL, cex=1, maxdepth=500, maxSNPs=1e5, res=300, color=NULL, ind=FALSE, ncores=1, ...){
                  config <- private$config[[1]]
                  if(any(is.na(config))) config[which(is.na(config))] <- private$config_infer[[1]][which(is.na(config))]
                  freq <- sapply(config, function(x) {
                    if(x == 1) return(c(0.25,0.5,0.25))
                    else if(x == 2 | x == 4) return(c(0,0.5,0.5))
                    else if(x == 3 | x == 5) return(c(0.5,0.5,0))
                  })
                  GUSbase::cometPlot(ref=private$ref[[1]], alt=private$alt[[1]], ploid=2, gfreq=freq, file=filename, cex=cex, 
                                     maxdepth=maxdepth, maxSNPs=maxSNPs, res=res, ind=ind, ncores = ncores, indID=private$indID, color=color, ...)
                  return(invisible())
                },
                rocketPlot = function(ploid=2, filename=NULL, cex=1, maxdepth=500, maxSNPs=1e5, res=300, scaled=TRUE, ...){
                  config <- private$config[[1]]
                  if(any(is.na(config))) config[which(is.na(config))] <- private$config_infer[[1]][which(is.na(config))]
                  freq <- sapply(config, function(x) {
                    if(x == 1) return(c(0.25,0.5,0.25))
                    else if(x == 2 | x == 4) return(c(0,0.5,0.5))
                    else if(x == 3 | x == 5) return(c(0.5,0.5,0))
                  })
                  GUSbase::rocketPlot(private$ref[[1]], private$alt[[1]], ploid=2, file=filename, gfreq=freq, cex=cex, 
                                      maxdepth=maxdepth, maxSNPs=maxSNPs, res=res, scaled=scaled, ...)
                  return(invisible())
                },
                # Ratio of alleles for heterozygous genotype calls (observed vs expected)
                RDDPlot = function(filename=NULL, maxdepth=500, maxSNPs=1e5, ...){
                  config <- private$config[[1]]
                  if(any(is.na(config))) config[which(is.na(config))] <- private$config_infer[[1]][which(is.na(config))]
                  freq <- sapply(config, function(x) {
                    if(x == 1) return(c(0.25,0.5,0.25))
                    else if(x == 2 | x == 4) return(c(0,0.5,0.5))
                    else if(x == 3 | x == 5) return(c(0.5,0.5,0))
                  })
                  GUSbase::RDDPlot(ref=private$ref[[1]], alt=private$alt[[1]], ploid=2, gfreq=freq, file=filename, maxdepth=maxdepth, maxSNPs=maxSNPs, ...)
                  return(invisible())
                },
                #### Write output
                writeLM = function(file, direct = "./", LG = NULL, what = NULL, inferGeno = TRUE){
                  ## do some checks
                  if(private$noFam == 1){
                    if(!is.vector(file) || !is.character(file) || length(file) != 1)
                      stop("Filename input is invalid")
                    if(!is.character(direct) || !is.vector(direct) || length(direct) != 1)
                      stop("Invalid input for the path to the directory where file is to be written.")
                    filename <- paste0(tail(strsplit(file,split=.Platform$file.sep)[[1]],1),"_GUSMap.txt")
                    outpath <- dts(normalizePath(direct, winslash=.Platform$file.sep, mustWork=T))
                    outfile = file.path(outpath,filename)
                    if(!file.create(outfile,showWarnings = F))
                      stop("Unable to create output file.")
                    ## check for NULL what input
                    if(is.null(what)){
                      if(!is.null(private$LG_map)) what = "map"
                      else if(!is.null(private$LG_mat_bi) || !is.null(private$LG_pat_bi)) what = "LG-bi"
                    }
                    ## Find out which LGs to write to file
                    if(what == "map"){
                      if(is.null(private$LG_map) && is.null(private$para))
                        stop("No linkage maps available to write to file.")
                      if(is.null(LG)){
                        LGlist <- private$LG_map
                        LG <- 1:length(private$LG_map)
                      } else if(isValue(LG, type="pos_integer", minv=1, maxv=length(private$LG_map)))
                        stop("LGs to write to file is invalid.")
                      else
                        LGlist <- private$LG_map[LG]
                      if(!is.vector(inferGeno) || length(inferGeno) != 1 || !is.logical(inferGeno) || is.na(inferGeno))
                        stop("Argument `inferGeno` is invalid. Should be a logical value")
                      if(length(unlist(LGlist)) == 0)
                        stop("No maps have been computed for the specified linkage groups.")
                      ## Infer genotypes if required
                      if(inferGeno){
                        nInd = unlist(private$nInd)
                        OPGPs = private$para$OPGP[LG]
                        #matgroups <- which(unlist(lapply(OPGPs[LG], function(x) all(x %in% c(1:4,9:12)))))
                        #patgroups <- which(unlist(lapply(OPGPs[LG], function(x) all(x %in% c(1:8)))))
                        geno = list()
                        for(lg in 1:length(LG)){
                          snps = LGlist[[lg]]
                          if(length(snps) == 0)
                            next
                          else{
                            nSnps = length(snps)
                            ref = private$ref[[1]][,snps]
                            alt = private$alt[[1]][,snps]
                            # output from viterbi_fs_err: 0 = 00, 1 = 10, 2 = 01, 3 = 11 (mat/pat)
                            ms = viterbi_fs_err(private$para$rf[[lg]], 
                                                rep(private$para$ep[[lg]], length.out=nSnps),
                                                nInd, nSnps, ref, alt,
                                                OPGPs[[lg]])
                            ## 1 = on paternal chrosome, 0 = maternal chromsome
                            infer_pat = (ms > 1) + 1
                            infer_mat = apply(ms, 2, function(x) x %in% c(1,3)) + 1
                            parHap = OPGPtoParHap(OPGPs[[lg]])
                            geno_temp = sapply(1:nInd, function(x) rbind(parHap[cbind(infer_mat[x,]+2,1:nSnps)], parHap[cbind(infer_pat[x,],1:nSnps)]), simplify = FALSE)
                            geno[[lg]] = t(do.call("rbind", geno_temp))
                          }
                        }
                      }
                      ## compute the information to go in file
                      nLG = length(LGlist)
                      LGno <- rep(LG, unlist(lapply(LGlist, length)))
                      LGindx <- sequence(lapply(LGlist, unlist(length)))
                      chrom <- private$chrom[unlist(LGlist)]
                      pos <- private$pos[unlist(LGlist)]
                      depth <- colMeans(private$ref[[1]][,unlist(LGlist)] + private$alt[[1]][,unlist(LGlist)])
                      config_temp <- private$config[[1]]
                      config_temp[!is.na(private$config_infer[[1]])] <- private$config_infer[[1]][!is.na(private$config_infer[[1]])] 
                      segType <- (c("BI","PI","PI","MI","MI","U"))[config_temp[unlist(LGlist)]]
                      temp <- private$ref[[1]][,unlist(LGlist)] > 0 | private$alt[[1]][,unlist(LGlist)] > 0
                      callrate <- colMeans(temp)
                      matgroups <- which(unlist(lapply(private$para$OPGP[LG], function(x) all(x %in% c(1:4,9:12)))))
                      patgroups <- which(unlist(lapply(private$para$OPGP[LG], function(x) all(x %in% c(1:8)))))
                      rf_m <- c(unlist(lapply(private$para$rf[matgroups], function(y) {if(!is.na(y[1])) return(format(round(cumsum(c(0,y)),6),digits=6,scientific=F))} )),
                                rep(0,length(unlist(unlist(LGlist[patgroups])))))
                      rf_p <- c(rep(0,length(unlist(unlist(LGlist[matgroups])))),
                                unlist(lapply(private$para$rf[patgroups], function(y) {if(!is.na(y[1])) return(format(round(cumsum(c(0,y)),6),digits=6,scientific=F))} )))
                      ep <- format(round(unlist(sapply(1:length(LG), function(x) rep(private$para$ep[[x]], length.out=length(LGlist[[x]])))),8),digits=8,scientific=F)
                      ## add geno information if required
                      if(inferGeno){
                        temp = do.call("rbind", geno)
                        temp2 = split(temp, rep(1:ncol(temp), each = nrow(temp)))
                        names(temp2) = paste0(c("MAT_","PAT_"), rep(private$indID[[1]], rep(2,nInd)))
                        out <- c(list(LG=LGno, LG_POS=LGindx, CHROM=chrom, POS=pos, TYPE=segType, RF_PAT=rf_p, RF_MAT=rf_m,
                                      ERR=ep, MEAN_DEPTH=depth, CALLRATE=callrate), temp2)
                      } else 
                        out <- list(LG=LGno, LG_POS=LGindx, CHROM=chrom, POS=pos, TYPE=segType, RF_PAT=rf_p, RF_MAT=rf_m,
                                    ERR=ep, MEAN_DEPTH=depth, CALLRATE=callrate)
                      ## write out file
                      data.table::fwrite(out, file = outfile, sep="\t", nThread = 1)  
                      cat(paste0("Linkage analysis results written to file:\nFilename:\t",filename,"\n")) 
                      return(invisible(NULL))
                    } else if(what == "LG-bi"){
                      if(is.null(private$LG_mat_bi) && is.null(private$LG_pat_bi))
                        stop("No ordered linkage groups available to write to file.")
                      LGlist <- c(private$LG_mat_bi,private$LG_pat_bi)
                      if(is.null(LG)){
                        LG <- 1:length(LGlist)
                      } else if(isValue(LG, type="pos_integer", minv=1, maxv=length(LGlist)))
                        stop("LGs to write to file is invalid.")
                      else
                        LGlist <- LGlist[LG]
                      ## compute the information to go in file
                      nLG = length(LGlist)
                      LGno <- rep(LG, unlist(lapply(LGlist, length)))
                      LGindx <- sequence(lapply(LGlist, unlist(length)))
                      chrom <- private$chrom[unlist(LGlist)]
                      pos <- private$pos[unlist(LGlist)]
                      depth <- colMeans(private$ref[[1]][,unlist(LGlist)] + private$alt[[1]][,unlist(LGlist)])
                      segType <- (c("BI","PI","PI","MI","MI","U"))[private$config[[1]][unlist(LGlist)]]
                      temp <- private$ref[[1]][,unlist(LGlist)] > 0 | private$alt[[1]][,unlist(LGlist)] > 0
                      callrate <- colMeans(temp)
                      out <- list(LG=LGno, LGPOS=LGindx, CHROM=chrom, POS=pos, TYPE=segType, MEAN_DEPTH=depth, CALLRATE=callrate)
                      ## write out file
                      data.table::fwrite(out, file = outfile, sep="\t", nThread = 1)  
                      cat(paste0("Linkage analysis results written to file:\nFilename:\t",filename,"\n")) 
                      return(invisible(NULL))
                    } else
                      stop("Invalid what input.")
                  }
                }
                ##############################################################
              ),
              private = list(
                config_orig  = NULL,
                config_infer_orig = NULL,
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
                LG_map       = NULL,
                LG_mat       = NULL,
                LG_mat_bi    = NULL,
                LG_pat       = NULL,
                LG_pat_bi    = NULL,
                summaryInfo  = NULL,
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
                        LODvalue = numeric(max(nLG,2))
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
              ), lock_objects = FALSE
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

#### Some functions from the kutils package for removing trailing spaces for filenames.
dts <- function (name)
  gsub("/$", "", dms(name))
dms <- function(name)
  gsub("(/)\\1+", "/", name)
