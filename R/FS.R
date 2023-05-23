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
#' FS object
#' 
#' Class for storing RA data and associated functions for analysis of full-sib family populations.
#' 
#' @usage
#' ## Create FS object
#' FSobj <- makeFS(RAobj, pedfile, family=NULL, inferSNPs = FALSE,
#'                 filter=list(MAF=0.05, MISS=0.2, BIN=100, DEPTH=5, PVALUE=0.01, MAXDEPTH=500))
#'
#' ## Functions (Methods) of an FS object
#' FSobj$addBIsnps(LODthres = 10, nComp = 10)
#' FSobj$addSNPs(LODthres = 10, nComp = 10)
#' FSobj$computeMap(chrom = NULL, init_r = 0.01, ep = 0.001, method = "optim", sexSpec = F, err = T, mapped = T, nThreads = 1)
#' FSobj$createLG(parent = "both", LODthres = 10, nComp = 10)
#' FSobj$maskSNP(snps)
#' FSobj$mergeLG(LG, where = NULL, mergeTo = NULL)
#' FSobj$orderLG(chrom = NULL, mapfun = "haldane", weight = "LOD2", ndim = 30, spar = NULL)
#' FSobj$plotChr(parent = "maternal", mat = "rf", filename = NULL, chrS = 2, lmai = 2)
#' FSobj$plotLG(parent = "maternal", LG = NULL, mat = "rf", filename = NULL, interactive = TRUE, what = NULL)
#' FSobj$plotLM(LG = NULL, fun = "haldane", col = "black")
#' FSobj$plotSyn()
#' FSobj$print(what = NULL, ...)
#' FSobj$removeLG(LG, where = NULL)
#' FSobj$removeSNP(snps, where = NULL)
#' FSobj$rf_2pt(nClust = 2, err = FALSE)
#' FSobj$unmaskSNP(snps)
#' FSobj$writeLM(file, direct = "./", LG = NULL, what = NULL, inferGeno=TRUE)
#' 
#' @details
#' An FS object is created from the \code{\link{makeFS}} function and contains RA data,
#' various statistics of the dataset that have been computed, and functions (methods)
#' for analyzing the data. Information in an FS object are specific to full-sib family populations.
#' 
#' @section Methods(Functions):
#' \describe{
#' \item{\code{\link{$addBIsnps}}}{Add Both-Informative (BI) snps to the maternal and/or
#' paternal linkage groups and merge the parental linkage groups.}
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
#' \item{\code{\link{$print}}}{Print summary information for FS object.}
#' \item{\code{\link{$removeLG}}}{Remove linkage group(s).}
#' \item{\code{\link{$removeSNP}}}{Remove SNP(s) from linkage group(s).}
#' \item{\code{\link{$rf_2pt}}}{Compute the 2-point recombination fraction (and LOD score) between all SNP pairs.}
#' \item{\code{\link{$unmaskSNP}}}{Unmask SNP(s) in the data set.}
#' \item{\code{\link{$writeLM}}}{Write linkage mapping results to a file.}
#' }
#' @format NULL
#' @author Timothy P. Bilton
#' @name FS
#' @export

### R6 class for creating a data format for full-sib families
FS <- R6::R6Class("FS",
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
                    if(!is.null(private$LG_mat) | !is.null(private$LG_pat))
                      what <- c(what, "LG-pts")
                    if(!is.null(private$LG))
                      what <- c(what, "LG-comb")
                    if(!is.null(private$para))
                      what <- c(what, "map")
                  } else if(!is.vector(what) || !is.character(what) ||
                            any(!(what %in% c("data","LG-pts","LG-comb","map"))))
                    stop("Argument for what output to print is invalid")
                  ## printe the required output
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
                    if(any(what == "LG-comb")){
                      if(is.null(private$LG))
                        warning("Linkage groups do not have any BI SNPs. Use the '$addBIsnps' to add BIsnps to the linkage groups.")
                      else{
                        cat("Combined Linkage Group Summary:\n")
                        MI <- unlist(lapply(private$LG, function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) %in% c(4,5), na.rm=T)))
                        PI <- unlist(lapply(private$LG, function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) %in% c(2,3), na.rm=T)))
                        BI <- unlist(lapply(private$LG, function(x) sum(c(private$config[[1]][x],private$config_infer[[1]][x]) == 1, na.rm=T)))
                        TOTAL <- unlist(lapply(private$LG, length))
                        tab <- cbind(LG=1:(length(private$LG)),MI,PI,BI,TOTAL)
                        tab <- stats::addmargins(tab, margin = 1)
                        rownames(tab) <- NULL
                        tab[nrow(tab), 1] <- "TOTAL"
                        prmatrix(tab, rowlab = rep("",nrow(tab)), quote=F)
                      }
                    }
                    if(any(what == "map")){
                      if(is.null(private$para))
                        warning("no maps have been estimated. Use the '$computeMap' function to compute some maps.")
                      else{
                        cat(private$summaryInfo$map[[1]])
                        #cat(private$summaryInfo$map[[2]])
                        prmatrix(private$summaryInfo$map[[2]], rowlab = rep("",nrow(private$summaryInfo$map[[2]])), quote=F)
                        #if(length(private$summaryInfo$map) == 5){
                        #  cat(private$summaryInfo$map[[4]])
                        #  prmatrix(private$summaryInfo$map[[5]], rowlab = rep("",nrow(private$summaryInfo$map[[5]])), quote=F)
                        #}
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
                    if(is.null(private$LG)) where = "LG-pts"
                    else where = "LG-comb"
                  }
                  if(where == "LG-pts"){
                    if (is.null(private$LG_mat) & is.null(private$LG_pat))
                      stop("Linakge groups have not been formed. Use the '$createLG' function to create linkage groups.")
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
                  else if(where == "LG-comb"){
                    if(is.null(private$LG))
                      stop("There are no combined linkage groups with BI SNPs. Use the '$addBIsnps' to create combined linkage groups.")
                    ## remove SNPs from the combined groups
                    for(lg in 1:length(private$LG)){
                      lgremove <- which(private$LG[[lg]] %in% snps)
                      if(length(lgremove) > 0)
                        private$LG[[lg]] <- private$LG[[lg]][-lgremove]
                    }
                    ## check for empty LGs
                    empty <- lapply(private$LG, length) == 0
                    if(any(empty))
                      private$LG[which(empty)] <- NULL
                    ## check whether all the linkage groups have been removed
                    if(length(private$LG) == 0){
                      message("All combined linkage groups removed. Use the `$addBIsnps` function to recreate the combined linkage groups.")
                      private$LG <- NULL
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
                    if(is.null(private$LG)) where = "LG-pts"
                    else where = "LG-comb"
                  }
                  if(where == "LG-comb"){
                    if(is.null(private$LG))
                      stop("There are no combined linkage groups with BI SNPs. Use the '$addBIsnps' to create combined linkage groups.")
                    nLG <- length(private$LG)
                    if(isValue(LG, type="pos_integer", minv=1,maxv=nLG))
                      stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nLG))
                    else{
                      private$LG <- private$LG[which(!(1:nLG %in% LG))]
                    }
                    ## check whether all the linkage groups have been removed
                    if(length(private$LG) == 0){
                      message("All combined linkage groups removed. Use the `$addBIsnps` function to recreate the combined linkage groups.")
                      private$LG <- NULL
                    }
                  }  else if(where == "LG-pts"){
                    if(is.null(private$LG_mat) || is.null(private$LG_pat))
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
                    if(is.null(private$LG)) where = "LG-pts"
                    else where = "LG-comb"
                  }
                  if(where == "LG-pts"){
                    ## Check that we have linkage groups
                    if(is.null(private$LG_mat) || is.null(private$LG_pat))
                      stop("No linkage groups exist. Please use the $createLG function to create some linkage groups")
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
                        private$LG_mat[[matgroups[1]]] <- mergedLG
                        private$LG_mat[matgroups[-1]] <- NULL
                        ## remove paternal linkage groups and change config
                        added_pat <- private$LG_pat[patgroups][[1]][private$LG_pat[patgroups][[1]] %in% c(private$group$MI,private$group$PI)]
                        if(length(added_pat) > 0)
                          private$config[[1]][added_pat] <- c(private$config[[1]][added_pat] %% 2) + 4
                        added_pat_infer <- private$LG_pat[patgroups][[1]][private$LG_pat[patgroups][[1]] %in% private$group_infer$SI]
                        if(length(added_pat_infer) > 0)
                          private$config_infer[[1]][added_pat_infer] <- c(private$config_infer[[1]][added_pat_infer] %% 2) + 4
                        private$LG_pat[patgroups] <- NULL
                      } else if(mergeTo == "paternal"){
                        private$LG_pat[[patgroups[1]]] <- mergedLG
                        private$LG_pat[patgroups[-1]] <- NULL
                        ## remove maternal linkage groups and change config
                        added_mat <- private$LG_mat[matgroups][[1]][private$LG_mat[matgroups][[1]] %in% c(private$group$MI,private$group$PI)]
                        if(length(added_mat) > 0)
                          private$config[[1]][added_mat] <- c(private$config[[1]][added_mat] %% 2) + 2
                        added_mat_infer <- private$LG_mat[matgroups][[1]][private$LG_mat[matgroups][[1]] %in% private$group_infer$SI]
                        if(length(added_mat_infer) > 0)
                          private$config_infer[[1]][added_mat_infer] <- c(private$config_infer[[1]][added_mat_infer] %% 2) + 2
                        private$LG_mat[matgroups] <- NULL
                      }
                    } else{
                      ## remove the LGs that were merged
                      if(length(matgroups) > 1)
                        private$LG_mat[matgroups[-1]] <- NULL
                      if(length(patgroups) > 1)
                        private$LG_pat[patgroups[-1]] <- NULL
                    }
                  } else if(where == "LG_comb"){
                    if(is.null(private$LG))
                      stop("There are no combined linkage groups with BI SNPs. Use the '$addBIsnps' to create combined linkage groups.")
                    nLG <- length(private$LG)
                    ## check input
                    if(isValue(LG, type="pos_integer", minv=1,maxv=nLG))
                      stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nLG))
                    LG <- unique(LG)
                    if(length(LG) == 1)
                      stop("Only one unique linkage group supplied. Need at least two linkage groups to merge")
                    ## merge LGs into the first
                    private$LG[[LG[1]]] <- unlist(private$LG[LG])
                    ## drop remaining groups
                    private$LG[LG[-1]] <- NULL
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
                    stop("Recombination fractions and LOD scores have not been computed.\nUse $rf_2pt() to compute the recombination fractions and LOD scores.")
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
                  
                  ## Reset the groups 
                  if(reset){
                    private$config <- private$config_orig
                    private$config_infer <- private$config_infer_orig
                    private$LG_mat = private$LG_pat = private$LG_mat_bi = private$LG_pat_bi = NULL
                    private$summaryInfo <- private$summaryInfo$data
                  } else if(parent == "maternal"){
                    private$LG_mat = private$LG_mat_bi = NULL
                  } else if(parent == "paternal"){
                    private$LG_pat = private$LG_pat_bi = NULL
                  } else if(parent == "both"){
                    private$LG_mat = private$LG_pat = private$LG_mat_bi = private$LG_pat_bi = NULL
                  }
                  
                  private$config <- private$config_orig
                  private$config_infer <- private$config_infer_orig
                  
                  
                  ## Create the groups
                  if(parent == "maternal" || parent == "both")
                    private$LG_mat <- createLG(private$group, private$LOD, "maternal", LODthres, nComp, masked=private$masked)
                  if(parent == "paternal" || parent == "both")
                    private$LG_pat <- createLG(private$group, private$LOD, "paternal", LODthres, nComp, masked=private$masked)
                  if(!is.null(private$LG))
                    private$LG <- NULL
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
                addBIsnps = function(LODthres=10, nComp=10, binsize=0){
                  
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
                  if(is.null(binsize) || GUSbase::checkVector(binsize, minv=0, type="pos_numeric", equal=FALSE) || length(binsize) != 1){
                    warning("Minimum distance between adjacent SNPs is not specified or is invalid. Setting to 0:")
                    binsize <- 0 
                  }
                  
                  # ## Use the bin information to work out which SNPs are on the same Tag and merge LGs this way
                  # if(binsize > 0){
                  #   ## determine bins
                  #   temp_bins = sapply(unique(private$chrom), function(x){
                  #     if(sum(private$chrom == x) == 1)
                  #       return(1)
                  #     else
                  #       return(cumsum(c(1,diff(private$pos[which(private$chrom == x)]) < binsize)))
                  #   })
                  #   bin_size = cumsum(unlist(lapply(temp_bins, function(y) length(unique(y)))))
                  #   bins = unlist(c(temp_bins[1], sapply(2:length(temp_bins), function(x){
                  #     temp_bins[[x]] + bin_size[x-1]
                  #   })))
                  #   
                  #   ## combine the LGs together:
                  #   for(lgmat in 1:length(private$LG_mat)){
                  #     snps = private$LG_mat[[lgmat]]
                  #     whichbins = bins[snps]
                  #     which
                  #     
                  #   }
                  #   
                  # }

                  ## Find the unmapped loci
                  unmapped <- sort(unlist(c(private$group$BI,private$group_infer$BI)))
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
                    if(length(newLGlist_mat) == 0 || length(newLGlist_pat) == 0){
                      cat("No BI SNPs mapped to maternal or/and paternal linkage groups. Un able to merge Linkage groups")
                    }
                    else{
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
                      indmat <- 1:nmat
                      indpat <- 1:npat
                      while(!finish){
                        if(length(mappedLG_mat) == nmat || length(mappedLG_pat) == npat){
                          finish = TRUE
                          next
                        }
                        indrow <- which(!(indmat %in% mappedLG_mat))
                        indcol <- which(!(indpat %in% mappedLG_pat))
                        tempBI <- matrix(tabBI[indrow,indcol], nrow=length(indrow), ncol=length(indcol))
                        maxv <- max(tempBI, na.rm=T)
                        if(maxv == 0){
                          finish = TRUE
                          next
                        }
                        temp <- which(tempBI == maxv, arr.ind=T)[1,]
                        tabBI[indrow[temp[1]], indcol[temp[2]]] <- NA
                        mappedLG_mat <- c(mappedLG_mat,indrow[temp[1]])
                        mappedLG_pat <- c(mappedLG_pat,indcol[temp[2]])
                      }
                      nmapped <- length(mappedLG_mat)
                      if(nmapped == 0)
                        stop("No BI SNPs could be mapped. Maybe try a different LOD threshold.")
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
                            LGmerged[[matlg]] <- c(LGmerged[[matlg]],newLGlist_pat[[mappedLG_pat[lg]]])
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
                  }
                  else if(exists("newLGlist_mat")){
                    private$LG <- newLGlist_mat
                  }
                  else if(exists("newLGlist_pat")){
                    private$LG <- newLGlist_pat
                  }
                  private$para <- NULL            ## reset the computed maps
                  private$summaryInfo$map <- NULL ## reset the computed maps summary
                  self$print(what = "LG-comb")
                  return(invisible(NULL))
                },
                ## Function for ordering linkage groups
                orderLG = function(LG = NULL, mapfun = "haldane", weight="LOD2", binsize=0, ndim=30, spar=NULL, filename=NULL){
                  ## do some checks
                  if(is.null(private$LG))
                    stop("There are no combined linkage groups. Please use the $addBISNPs to create combined LGs")
                  if(is.null(LG))
                    LG <- 1:length(private$LG)
                  else if(!is.vector(LG) && !is.numeric(LG) && any(LG < 0) && LG > length(private$LG))
                    stop("Invalid Linkage Group input")
                  if(!(mapfun %in% c("morgan","haldane","kosambi")))
                    stop("Unknown mapping function")
                  if(!(weight %in% c("LOD","LOD2","none")))
                    stop("Unknown weighting function")
                  if(is.null(binsize) || GUSbase::checkVector(binsize, minv=0, type="pos_numeric", equal=FALSE) || length(binsize) != 1){
                    warning("Minimum distance between adjacent SNPs is not specified or is invalid. Setting to 0:")
                    binsize <- 0 
                  }
                  if(!is.null(filename) && (!is.vector(filename) || !is.character(filename) || length(filename) != 1))
                    stop("Specified filename is invalid")
                  if(!is.null(filename)){
                    temp <- file.create(paste0(filename,"_LG",LG[1],".pdf"))
                    if(!temp)
                      stop("Error in filename. Cannot create pdf")
                  }
                  nChr = length(LG)
                  tmp_config = private$config[[1]]
                  tmp_config[which(is.na(tmp_config))] = private$config_infer[[1]][!is.na(private$config_infer[[1]])]
                  ## if binning SNPs then work out the bins
                  if(binsize > 0){
                    ## determine bins
                    temp_bins = sapply(unique(private$chrom), function(x){
                      if(sum(private$chrom == x) == 1)
                        return(1)
                      else
                        return(cumsum(c(1,diff(private$pos[which(private$chrom == x)]) > binsize)))
                    })
                    bin_size = cumsum(unlist(lapply(temp_bins, function(y) length(unique(y)))))
                    bins = unlist(c(temp_bins[1], sapply(2:length(temp_bins), function(x){
                      temp_bins[[x]] + bin_size[x-1]
                    })))
                    tmp_depth = private$ref[[1]] + private$alt[[1]]
                  }
                  ## order the linkage groups
                  for(lg in LG){
                    ind <- private$LG[[lg]]
                    ind_mi <- which(ind %in% which(tmp_config %in% c(4,5)))
                    ind_pi <- which(ind %in% which(tmp_config %in% c(2,3)))
                    ind_bi = which(ind %in% which(tmp_config %in% 1))
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
                    ## if we are combining information from SNPs on the same tag
                    if(binsize > 0){
                      bins_mi = bins[ind[ind_mi]]
                      bins_pi = bins[ind[ind_pi]]
                      bin_index = bins[ind]
                      #all_bins = unique(bin_index)
                      bins_in_common = names(which(table(bin_index)>1))
                      #bins_in_common = all_bins[which((all_bins %in% bins_mi) & (all_bins %in% bins_pi))]
                      newmat.indx = which(!(bin_index %in% bins_in_common))
                      ## now iterate over the bins to combine
                      for(bin in bins_in_common){
                        tmp2 = which(bin_index == bin)
                        tmp = ind[tmp2]
                        tag_pi = tmp_config[tmp] %in% c(2,3)
                        tag_mi = tmp_config[tmp] %in% c(4,5)
                        tag_bi = tmp_config[tmp] %in% 1
                        if(any(tag_mi) & any(tag_pi)){
                          tmp_indx_mi = which(tag_mi)[which.max(colMeans(tmp_depth[,tmp][,which(tag_mi), drop=FALSE]))]
                          tmp_indx_pi = which(tag_pi)[which.max(colMeans(tmp_depth[,tmp][,which(tag_pi), drop=FALSE]))]
                          wmat[tmp2[tmp_indx_mi],ind_pi] = wmat[tmp2[tmp_indx_pi],ind_pi]
                          wmat[tmp2[tmp_indx_mi],ind_bi] = colMeans(wmat[tmp2[c(tmp_indx_mi, tmp_indx_pi)],ind_bi])
                          dmat[tmp2[tmp_indx_mi],ind_pi] = dmat[tmp2[tmp_indx_pi],ind_pi]
                          dmat[tmp2[tmp_indx_mi],ind_bi] = colMeans(dmat[tmp2[c(tmp_indx_mi, tmp_indx_pi)],ind_bi])
                          newmat.indx = c(newmat.indx, tmp2[tmp_indx_mi])
                        } else if(any(tag_bi)){
                          tmp_indx_bi = which(tag_bi)[which.max(colMeans(tmp_depth[,tmp][,which(tag_bi), drop=FALSE]))]
                          newmat.indx = c(newmat.indx, tmp2[tmp_indx_bi])
                        } else if(all(tag_mi)){
                          newmat.indx = c(newmat.indx, tmp2[which.max(colMeans(tmp_depth[,tmp][,which(tag_mi), drop=FALSE]))])
                        } else if(all(tag_pi)){
                          newmat.indx = c(newmat.indx, tmp2[which.max(colMeans(tmp_depth[,tmp][,which(tag_pi), drop=FALSE]))])
                        } else if(all(tag_bi)){
                          newmat.indx = c(newmat.indx, tmp2[which.max(colMeans(tmp_depth[,tmp][,which(tag_bi), drop=FALSE]))])
                        } else
                          stop("Error in code.")
                      }
                      tmp.wmat = wmat[newmat.indx,newmat.indx]
                      tmp.dmat = dmat[newmat.indx,newmat.indx]
                      tmp.binID = bins[ind[newmat.indx]]
                      ## unconstrainted MDS
                      MDS <- smacof::smacofSym(delta=tmp.dmat, ndim=ndim, weightmat = tmp.wmat, itmax = 1e+05)
                      ## principal curve
                      pcurve <- princurve::principal_curve(MDS$conf, maxit = 150, spar = spar)
                      neworder = pcurve$ord
                      neworder = order(as.numeric(factor(bin_index, levels = tmp.binID[neworder])))
                    } else{ ## original approach
                      ## unconstrainted MDS
                      MDS <- smacof::smacofSym(delta=dmat, ndim=ndim, weightmat = wmat, itmax = 1e+05)
                      ## principal curve
                      pcurve <- princurve::principal_curve(MDS$conf, maxit = 150, spar = spar)
                      neworder = pcurve$ord
                    }
                    ## plot the results
                    if(is.null(filename)) temp_par <- par(no.readonly=T)
                    else grDevices::pdf(paste0(filename,"_LG",lg,".pdf"), width=8, height=8)
                    graphics::par(mfrow = c(1,2))
                    graphics::plot(MDS$conf[,1],MDS$conf[,2], ylab="Dimension 2", xlab="Dimension 1", type="n")
                    graphics::text(MDS$conf, labels=ind)
                    graphics::lines(pcurve)
                    graphics::image(private$rf[ind,ind][neworder,neworder], axes=F,
                                    col=grDevices::heat.colors(100))
                    if(is.null(filename)) graphics::par(temp_par)
                    else grDevices::dev.off()
                    ## Set the new order
                    private$LG[[lg]] <- ind[neworder]
                  }
                  ## Acknowledge paper we are using algorithm from
                  cat("Using the MDS scaling algorithm to order SNPs.\n",
                      "Please cite:\n",
                      "\tPreedy & Hackett (2016). Theor Appl Genet. 129:2117-2132\n",sep="")
                  return(invisible(NULL))
                },
                ## Function for setting temp LGs to new LGs
                # setLG = function(parent="both"){
                #   if(parent == "both"){
                #     private$LG_mat <- private$LG_mat_temp
                #     private$LG_mat_temp <- NULL
                #     private$LG_pat <- private$LG_pat_temp
                #     private$LG_pat_temp <- NULL
                #     cat("Temporary linkage groups now set as new linkage groups\n")
                #   }
                #   else if(parent == "maternal"){
                #     private$LG_mat <- private$LG_mat_temp
                #     private$LG_mat_temp <- NULL
                #     cat("Temporary maternal linkage groups now set as new linkage groups\n")
                #   }
                #   else if(parent == "paternal"){
                #     private$LG_pat <- private$LG_pat_temp
                #     private$LG_pat_temp <- NULL
                #     cat("Temporary paternal linkage groups now set as new linkage groups\n")
                #   }
                #   return(invisible())
                # },
                ## Function for plotting linkage groups
                plotLG = function(parent = "maternal", LG=NULL, mat="rf", filename=NULL, interactive=TRUE, what = NULL, ...){
                  ## do some checks
                  if(!is.vector(mat) || !is.character(mat) || length(mat) != 1 || !(mat %in% c('rf','LOD')))
                    stop("Argument specifying which matrix to plot (argument 1) must be either 'rf' or 'LOD'")
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop("parent argument is not a string of length one or is invalid:
Please select one of the following:
  maternal: Add BI SNPs to MI LGs
  paternal: Add BI SNPs to PI LGs
  both:     Add BI SNPs to both MI and PI LGs")
                  if(!is.vector(interactive) || length(interactive) != 1 || !is.logical(interactive))
                    stop("Argument sepcifting whether to produce an interactive heatmap or not is invalid.")
                  if(!is.null(filename) && (!is.vector(filename) || !is.character(filename) || length(filename) != 1))
                    stop("Specified filename is invalid")
                  ## Check the what input
                  if(is.null(what)){
                    if(is.null(private$LG)) what = "LG-pts"
                    else what = "LG-comb"
                  }
                  plotly.arg <- list(...)
                  if(is.null(plotly.arg$selfcontained))
                    plotly.arg$selfcontained = FALSE
                  
                  if(private$noFam == 1){
                    ## Work out which LGs list to use
                    if(what == "LG-pts"){
                      if(parent == "both"){
                        if(is.null(private$LG_pat) || is.null(private$LG_mat))
                          stop("No linkage groups are available to be plotted. Please use the $createLG function to create some linkage groups")
                        else
                          LGlist <- c(private$LG_mat,private$LG_pat)
                      } else if(parent == "maternal"){
                        if(is.null(private$LG_mat))
                          stop("No maternal linkage groups are available to be plotted. Please use the $createLG function to create some maternal linkage groups")
                        else
                          LGlist <- private$LG_mat
                      } else{
                        if(is.null(private$LG_pat))
                          stop("No paternal linkage groups are available to be plotted. Please use the $createLG function to create some paternal linkage groups")
                        else
                          LGlist <- private$LG_pat
                      }
                    } else if(what == "LG-comb"){
                      if(is.null(private$LG))
                        stop("There are no combined linkage groups with BI SNPs. Use the '$addBIsnps' to create combined linkage groups.")
                      else
                        LGlist <- private$LG
                      ## Check which type of SNPs we are plotting
                      if(parent == "maternal"){
                        mi_ind <- lapply(LGlist, function(x) x[which(x %in% c(which(private$config[[1]] %in% c(1,4,5)),
                                                                              which(private$config_infer[[1]] %in% c(1,4,5))))])
                        LGlist <- mi_ind[which(unlist(lapply(mi_ind, length))!=0)]
                      } else if (parent == "paternal"){
                        pi_ind <- lapply(LGlist, function(x) x[which(x %in% c(which(private$config[[1]] %in% c(1:3)),
                                                                              which(private$config_infer[[1]] %in% c(1:3))))])
                        LGlist <- pi_ind[which(unlist(lapply(pi_ind, length))!=0)]
                      }
                    } else
                      stop("invalid argument 'what'.")
                    ## Work out if we want a subset of the LGs
                    if(!is.null(LG)){
                      if(isValue(LG, type="pos_integer", minv=1, maxv=length(LGlist)))
                        stop(paste0("At least one linkage group number does not exist. Indices must be between 1 and",length(LGlist),"\n"))
                      LGlist <- LGlist[LG]
                      names(LGlist) <- LG
                    }
                    else
                      names(LGlist) <- 1:length(LGlist)
                    

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
                      hovertext <- matrix(paste(matrix(paste0("row: ",chrom.ind), nrow=nn, ncol=nn), 
                                                matrix(paste0("col: ",chrom.ind), nrow=nn, ncol=nn, byrow=T), paste0("rf: ",round(temprf,4)), sep="<br>"),
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
                          grDevices::png(paste0(filename,".png"), width=nn+1,height=nn+1)
                      }
                      else
                        temp_par <- graphics::par(no.readonly = TRUE)
                      graphics::par(mar=rep(0,4),oma=c(0,0,0,0), mfrow=c(1,1), xaxt='n',yaxt='n',bty='n',ann=F)
                      graphics::image(temprf, x=1:nn, y=1:nn, zlim=zlim, col=mycol)
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
                plotChr = function(chrom = NULL, parent = "maternal", mat=c("rf"), filename=NULL, chrS=2, lmai=2){
                  ## do some checks
                  if(!is.vector(mat) || !is.character(mat) || length(mat) != 1 || !(mat %in% c('rf','LOD')))
                    stop("Argument specifying which matrix to plot (argument 1) must be either 'rf' or 'LOD'")
                  if(!is.character(parent) || length(parent) != 1 || !(parent %in% c("maternal","paternal","both")))
                    stop(paste("parent argument is not a string of length one or is invalid:",
                               "Please select one of the following:",
                               "maternal: Add BI SNPs to MI LGs",
                               "paternal: Add BI SNPs to PI LGs",
                               "both:     Add BI SNPs to both MI and PI LGs", sep="\n"))
                  
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
                      LG <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,4,5))), simplify=F)
                    else if(parent == "paternal")
                      LG <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,2,3))), simplify=F)
                    else if(parent == "both")
                      LG <- sapply(names, function(x) which((private$chrom == x) & !private$masked & (private$config[[1]] %in% c(1,2,3,4,5))), simplify=F)
                    ## plot the chromsomes rf info
                    if(mat == "rf")
                      plotLG(mat=private$rf, LG=LG[chrom], filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=T, type="rf")
                    else if(mat == "LOD")
                      plotLG(mat=private$LOD, LG=LG[chrom], filename=filename, names=names, chrS=chrS, lmai=lmai, chrom=T, type="LOD")
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
                  if(is.null(private$para))
                    stop("There are no linkage maps availabe to plot. Use the '$computeMap' function to compute linkage maps.")
                  nLGs = length(private$para[[1]])
                  if(is.null(LG))
                    LG <- 1:length(private$para$rf_p)
                  if(isValue(LG, type="pos_integer", minv=1, maxv=nLGs))
                    stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nLGs))
                  LG = LG[which(!unlist(lapply(private$para$rf_p[LG], function(x) any(is.na(x)))))]
                  nLGs = length(LG)
                  if(nLGs == 0)
                    stop("There are no linkage maps for the specified linkage groups.")
                  if(!is.vector(fun) || length(fun) != 1 || any(!(fun %in% c("morgan", "haldane", "kosambi"))))
                    stop("Input argument for mapping distance is invalid")
                  ## check the color input
                  if(!is.vector(is.character(col)) || !is.character(col))
                    stop("Input argument for colors of the linkage maps is invalid.")
                  else{
                    tc <- unname(sapply(col, function(y) tryCatch(is.matrix(col2rgb(y)), error = function(e) FALSE)))
                    if(any(!tc))
                      stop("At least one color in the color list is undefined")
                    else
                      col = rgb(t(col2rgb(col)), max=255)
                    if(length(col) == 1)
                      col <- rep(col,nLGs)
                    else if(length(col) != nLGs)
                      stop("Number of colors specified does not equal the number of linkage groups")
                  }
                  temp_par <- par(no.readonly = TRUE) # save the current plot margins
                  
                  ellipseEq_pos <- function(x) c2 + sqrt(b^2*(1-round((x-c1)^2/a^2,7)))
                  ellipseEq_neg <- function(x) c2 - sqrt(b^2*(1-round((x-c1)^2/a^2,7)))
                  mapDist = lapply(private$para$rf_p[LG],mfun, fun=fun, centi=TRUE)
                  yCoor = max(unlist(lapply(mapDist,sum)))
                  graphics::par(mar=rep(2,4), mfrow=c(1,1))
                  graphics::plot(NULL,NULL, ylim=c(yCoor+0.01*yCoor,-0.01*yCoor),xlim=c(0,0.25*(nLGs+1)), xaxt='n',bty='n',xlab="",ylab="")
                  err = 0.05
                  ## plot each map for each linkage group
                  for(lg in 1:nLGs){
                    cent = 0.25*lg
                    segments(cent-err,c(0,cumsum(mapDist[[lg]])),cent+err,c(0,cumsum(mapDist[[lg]])),col=col[lg])
                    segments(cent-err,-0.01*yCoor,cent-err,sum(mapDist[[lg]])+0.01*yCoor,col=col[lg])
                    segments(cent+err,-0.01*yCoor,cent+err,sum(mapDist[[lg]])+0.01*yCoor,col=col[lg])
                    ## Plot the curves at the ends
                    a=err; b=0.02*yCoor; c1=cent; c2=-0.01*yCoor;
                    curve(ellipseEq_neg, from=cent-err,to=cent+err,add=T,col=col[lg])
                    c2=sum(mapDist[[lg]])+0.01*yCoor;
                    curve(ellipseEq_pos, from=cent-err,to=cent+err,add=T,col=col[lg])
                  }
                  graphics::par(temp_par) # reset the plot margins
                  #if(!is.null(names))
                  #  mtext(names, side=1, at=seq(0.25,0.25*nLGs,0.25))
                },
                ## compare order with original and ordered markers
                plotSyn = function(){
                  if(is.null(private$LG))
                    stop("There are no combined linkage groups availabe. Use the '$addBIsnps' function to compute combined linkage groups.")
                  ## work out the LG and original orders
                  LGorder <- unlist(private$LG)
                  orgOrder <- sort(LGorder)
                  ## determine the breask on the plot
                  temp <- cumsum(unlist(lapply(private$LG, length)))
                  LGbreaks <- orgOrder[temp[-length(temp)]] + 0.5
                  chrBreaks <- which(diff(as.numeric(factor(private$chrom, levels=unique(private$chrom))))==1) + 0.5
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
                  graphics::mtext(text = 1:length(private$LG), side = 1, 
                                  at = apply(cbind(c(min(LGorder),LGbreaks),c(LGbreaks,max(LGorder))),1,mean))
                  graphics::par(temp_par) # reset the plot margins
                }, 
                ## Function for computing the rf's for each chromosome 
                computeMap = function(chrom=NULL, init_r=0.001, ep=0.001, method=NULL, sexSpec=FALSE, err=TRUE, multiErr=FALSE, 
                                      mapped=TRUE, nThreads=1, inferOPGP=TRUE, rfthres=0.1){
                  ## do some checks
                  if( !is.null(init_r) & !is.numeric(init_r) )
                    stop("Starting values for the recombination fraction needs to be a numeric vector or integer or a NULL object")
                  if( (length(ep) != 1 || !is.numeric(ep) || (ep <= 0 | ep >= 1)) )
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
                    list_chrom <- unique(private$chrom)
                    nchrom <- length(list_chrom)
                    if(!is.null(chrom) && !is.vector(chrom))
                      stop("chromosomes names must be a character vector.")
                    if(is.null(chrom)) # compute distances for all chromosomes
                      chrom <- list_chrom
                    else if(!(any(chrom %in% list_chrom)))
                      stop("At least one chromosome not found in the data set.")
                    else
                      chrom <- as.character(chrom) # ensure that the chromsome names are characters
                    wchrom <- which((list_chrom %in% chrom) & (list_chrom %in% unique(private$chrom[!private$masked])))
                    cat("Computing recombination fractions:\n")
                    if(is.null(private$para))
                      private$para <- list(OPGP=vector(mode = "list",length = nchrom),rf_p=vector(mode = "list",length = nchrom),
                                           rf_m=vector(mode = "list",length = nchrom),ep=vector(mode = "list",length = nchrom),
                                           loglik=vector(mode = "list",length = nchrom))
                    else if(length(private$para$rf_p) != nchrom)
                      private$para <- list(OPGP=vector(mode = "list",length = nchrom),rf_p=vector(mode = "list",length = nchrom),
                                           rf_m=vector(mode = "list",length = nchrom),ep=vector(mode = "list",length = nchrom),
                                           loglik=vector(mode = "list",length = nchrom))
                     for(i in wchrom){
                      cat("chromosome: ",list_chrom[i],"\n")
                      indx_chrom <- which((private$chrom == list_chrom[i]) & !private$masked)
                      ref_temp <- lapply(private$ref, function(x) x[,indx_chrom])
                      alt_temp <- lapply(private$alt, function(x) x[,indx_chrom])
                      ## estimate OPGP's
                      curOPGP <- private$para$OPGP[i]
                      if(inferOPGP || !is.list(curOPGP) || length(curOPGP[[1]]) != length(indx_chrom)){
                        tempOPGP <- list()
                        for(fam in 1:private$noFam){
                          tempOPGP <- c(tempOPGP,list(as.integer(infer_OPGP_FS(ref_temp[[fam]],alt_temp[[fam]],private$config[[fam]][indx_chrom], method="EM", nThreads=nThreads))))
                        }
                        private$para$OPGP[i] <- tempOPGP
                      }
                      ## estimate the rf's
                      if(!is.null(method))
                        MLE <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                         sexSpec=sexSpec, seqErr=err, method=method, nThreads=nThreads, multiErr=multiErr)
                      else{
                        MLE_init <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                              sexSpec=sexSpec, seqErr=err, method="EM", nThreads=nThreads, multiErr=multiErr, maxit=20)
                        MLE <- rf_est_FS(init_r=MLE_init$rf, ep=MLE_init$ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                         sexSpec=sexSpec, seqErr=err, method="optim", nThreads=nThreads, multiErr=multiErr)
                      }
                      
                      if(sexSpec){
                        private$para$rf_p[i]   <- list(MLE$rf_p)
                        private$para$rf_m[i]   <- list(MLE$rf_m)
                      } else{
                        private$para$rf_p[i]   <- list(MLE$rf)
                        private$para$rf_m[i]   <- list(MLE$rf)
                      }
                      private$para$ep[i]    <- list(list(MLE$ep))
                      private$para$loglik[i]<- list(MLE$loglik)
                      ### Check for SNPs with really large rf values
                      if(!sexSpec){
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
                  }
                  else{
                    ## Check the input
                    if((length(private$LG) == 0) && length(is.null(private$LG) == 0))
                      stop("No linkage groups have been formed.")
                    nchrom <- length(private$LG)
                    if(is.null(chrom)) # compute distances for all chromosomes
                      chrom <- 1:nchrom
                    else if(!is.vector(chrom) || chrom < 0 || chrom > length(private$LG) || round(chrom) != chrom)
                      stop("Linkage group input must be an integer vector between 1 and the number of linkage groups")
                    ## Compute rf's
                    cat("Computing recombination fractions:\n")
                    if(is.null(private$para))
                      private$para <- list(OPGP=replicate(nchrom, NA, simplify=FALSE),rf_p=replicate(nchrom, NA, simplify=FALSE),
                                           rf_m=replicate(nchrom, NA, simplify=FALSE),ep=replicate(nchrom, NA, simplify=FALSE),
                                           loglik=replicate(nchrom, NA, simplify=FALSE))
                    else if(length(private$para$rf_p) != length(private$LG))
                      private$para <- list(OPGP=replicate(nchrom, NA, simplify=FALSE),rf_p=replicate(nchrom, NA, simplify=FALSE),
                                           rf_m=replicate(nchrom, NA, simplify=FALSE),ep=replicate(nchrom, NA, simplify=FALSE),
                                           loglik=replicate(nchrom, NA, simplify=FALSE))
                    if(is.null(private$LG_map))
                      private$LG_map <- vector(mode="list", length = nchrom)
                    for(i in chrom){
                      cat("Linkage Group: ",i,"\n")
                      indx_chrom <- private$LG[[i]]
                      ref_temp <- lapply(private$ref, function(x) x[,indx_chrom])
                      alt_temp <- lapply(private$alt, function(x) x[,indx_chrom])
                      config_temp <- lapply(private$config, function(x) x[indx_chrom])
                      for(fam in 1:private$noFam){
                        tempIndx <- !is.na(private$config_infer[[fam]][indx_chrom])
                        config_temp[[fam]][tempIndx] <- private$config_infer[[fam]][indx_chrom[tempIndx]]
                      }
                      
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
                                         sexSpec=sexSpec, seqErr=err, method=method, nThreads=nThreads, multiErr=multiErr)
                      else{
                        MLE_init <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                              sexSpec=sexSpec, seqErr=err, method="EM", nThreads=nThreads, multiErr=multiErr, maxit=20)
                        MLE <- rf_est_FS(init_r=MLE_init$rf, ep=MLE_init$ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                         sexSpec=sexSpec, seqErr=err, method="optim", nThreads=nThreads, multiErr=multiErr)
                      }
                      if(sexSpec){
                        private$para$rf_p[i]   <- list(MLE$rf_p)
                        private$para$rf_m[i]   <- list(MLE$rf_m)
                      } else{
                        private$para$rf_p[i]   <- list(MLE$rf)
                        private$para$rf_m[i]   <- list(MLE$rf)
                      }
                      private$para$ep[i]       <- list(MLE$ep)
                      private$para$loglik[i]   <- list(MLE$loglik)
                      private$LG_map[[i]]      <- private$LG[[i]]
                      ### Check for SNPs with really large rf values
                      if(!sexSpec){
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
                    MI <- unlist(lapply(private$LG_map, function(x) sum(x %in% private$group$MI)))
                    PI <- unlist(lapply(private$LG_map, function(x) sum(x %in% private$group$PI)))
                    BI <- unlist(lapply(private$LG_map, function(x) sum(x %in% private$group$BI)))
                    TOTAL <- unlist(lapply(private$LG_map, length))
                    tab <- cbind(LG=1:nchrom,MI,PI,BI,TOTAL)
                    DIST_MAT <- extendVec(unlist(lapply(private$para$rf_m, function(x) round(sum(mfun(x, fun="haldane", centiM = T)),2))), nrow(tab))
                    DIST_PAT <- extendVec(unlist(lapply(private$para$rf_p, function(x) round(sum(mfun(x, fun="haldane", centiM = T)),2))), nrow(tab))
                    ERR <- extendVec(round(unlist(lapply(private$para$ep,mean)),5), nrow(tab))
                    tab <- cbind(tab,DIST_MAT,DIST_PAT,ERR)
                    tab <- stats::addmargins(tab, margin = 1, FUN=function(x) sum(x, na.rm=T))
                    rownames(tab) <- NULL
                    tab[nrow(tab), 1] <- "TOTAL"
                    if(ncol(tab) > 7)
                      tab[nrow(tab),8] <- ""
                    private$summaryInfo$map <- list("Linkage Map Summaries:\n", tab)
                    self$print(what = "map")
                    return(invisible(NULL))
                  }
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
                      else if(!is.null(private$LG)) what = "LG-comb"
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
                        geno = list()
                        for(lg in LG){
                          snps = LGlist[[lg]]
                          if(length(snps) == 0)
                            next
                          else{
                            nSnps = length(snps)
                            ref = private$ref[[1]][,snps]
                            alt = private$alt[[1]][,snps]
                            # output from viterbi_fs_err: 0 = 00, 1 = 10, 2 = 01, 3 = 11 (mat/pat)
                            ms = viterbi_fs_err_ss(private$para$rf_p[[lg]], private$para$rf_m[[lg]],
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
                      rf_p <- unlist(lapply(private$para$rf_p[LG], function(y) {if(!is.na(y[1])) return(format(round(cumsum(c(0,y)),6),digits=6,scientific=F))} ))
                      rf_m <- unlist(lapply(private$para$rf_m[LG], function(y) {if(!is.na(y[1])) return(format(round(cumsum(c(0,y)),6),digits=6,scientific=F))} ))
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
                      cat(paste0("Linkage analysis results written to file:\nFilename:\t",filename)) 
                      return(invisible(NULL))
                    } else if(what == "LG-comb"){
                      if(is.null(private$LG))
                        stop("No combined linkage groups available to write to file.")
                      if(is.null(LG)){
                        LGlist <- private$LG
                        LG <- 1:length(private$LG)
                      } else if(GUSbase::checkVector(LG, type="pos_integer", minv=1, maxv=length(private$LG)))
                        stop("LGs to write to file is invalid.")
                      else
                        LGlist <- private$LG[LG]
                      ## compute the information to go in file
                      nLG = length(LGlist)
                      LGno <- rep(LG, unlist(lapply(LGlist, length)))
                      LGindx <- sequence(lapply(LGlist, unlist(length)))
                      chrom <- private$chrom[unlist(LGlist)]
                      pos <- private$pos[unlist(LGlist)]
                      depth <- colMeans(private$ref[[1]][,unlist(LGlist)] + private$alt[[1]][,unlist(LGlist)])
                      segType <- (c("BI","PI","PI","MI","MI","U"))[private$config[[1]][unlist(private$LG[LG])]]
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
                },
                #### Diagonostic functions ####
                # Ratio of alleles for heterozygous genotype calls (observed vs expected)
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
                LG           = NULL,
                LG_map       = NULL,
                LG_mat       = NULL,
                LG_mat_temp  = NULL,
                LG_pat       = NULL,
                LG_pat_temp  = NULL,
                summaryInfo  = NULL,
                simPara      = NULL,
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
