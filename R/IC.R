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
#' IC object
#' 
#' Class for storing RA data and associated functions for analysis of intercross (F2) populations.
#' 
#' @usage
#' ## Create IC object
#' ICobj <- makeIC(RAobj, samID=NULL,
#'                 filter=list(MAF=0.05, MISS=0.2, BIN=100, DEPTH=5, PVALUE=0.01, MAXDEPTH=500))
#'
#' ## Functions (Methods) of an IC object
#' ICobj$computeMap(chrom = NULL, init_r = 0.01, ep = 0.001, method = "EM", err = T, mapped = T, nThreads = 1)
#' ICobj$createLG(LODthres = 10, nComp = 10)
#' ICobj$maskSNP(snps)
#' ICobj$mergeLG(LG)
#' ICobj$orderLG(chrom = NULL, mapfun = "haldane", weight = "LOD2", ndim = 30, spar = NULL)
#' ICobj$plotChr(mat = "rf", filename = NULL, chrS = 2, lmai = 2)
#' ICobj$plotLG(LG = NULL, mat = "rf", filename = NULL, interactive = TRUE)
#' ICobj$plotLM(LG = NULL, fun = "haldane", col = "black")
#' ICobj$plotSyn()
#' ICobj$print(what = NULL, ...)
#' ICobj$removeLG(LG)
#' ICobj$removeSNP(snps)
#' ICobj$rf_2pt(nClust = 2, err=FALSE)
#' ICobj$unmaskSNP(snps)
#' ICobj$writeLM(file, direct = "./", LG = NULL, what = NULL, inferGeno = TRUE)
#' 
#' @details
#' An IC object is created from the \code{\link{makeIC}} function and contains RA data,
#' various statistics of the dataset that have been computed, and functions (methods)
#' for analyzing the data. Information in an IC object are specific to intercross (F2) populations.
#' 
#' @section Methods(Functions):
#' \describe{
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
#' @name IC
#' @export

### R6 class for creating a data format for intercross (F2) population
IC <- R6::R6Class("IC",
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
                        if(!is.null(private$LG))
                          what <- c(what, "LG")
                        if(!is.null(private$para))
                          what <- c(what, "map")
                      } else if(!is.vector(what) || !is.character(what) ||
                                any(!(what %in% c("data","LG","map"))))
                        stop("Argument for what output to print is invalid")
                      ## printe the required output
                      if(private$noFam == 1){
                        if(any(what == "data"))
                          cat(private$summaryInfo$data, sep="")
                        if(any(what == "LG")){
                          if(is.null(private$LG))
                            warning("Linkage groups have not been formed. Use the '$createLG' function to form linkage groups.")
                          else{
                            cat("Combined Linkage Group Summary:\n")
                            BI <- unlist(lapply(private$LG, function(x) sum(x %in% private$group$BI)))
                            TOTAL <- unlist(lapply(private$LG, length))
                            tab <- cbind(LG=1:(length(private$LG)),BI,TOTAL)
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
                            prmatrix(private$summaryInfo$map[[2]], rowlab = rep("",nrow(private$summaryInfo$map[[2]])), quote=F)
                          }
                        }
                      }
                      else
                        cat("Not yet implemented")
                      return(invisible(NULL))
                    },
                    #############################################################
                    ## Function for removing SNPs from the linkage groups
                    removeSNP = function(snps){
                      ## some checks
                      if( !is.vector(snps) || !is.numeric(snps)  || any(is.na(snps)) || !all(snps == round(snps)) || any(snps < 1) || any(snps > private$nSnps) )
                        stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                      snps <- unique(snps) ## make sure no double ups in SNPs
                      if(is.null(private$LG))
                        stop("Linkage groups have not been formed. Use the '$createLG' function to create linkage groups.")
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
                        message("All linkage groups have been removed.")
                        private$LG <- NULL
                      }
                      return(invisible(NULL))
                    },
                    ## Function for removing linkage groups 
                    removeLG = function(LG){
                      if(is.null(private$LG))
                        stop("Linkage groups have not been formed. Use the '$createLG' function to create linkage groups.")
                      nLG <- length(private$LG)
                      if(isValue(LG, type="pos_integer", minv=1,maxv=nLG))
                        stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nLG))
                      else{
                        private$LG <- private$LG[which(!(1:nLG %in% LG))]
                      }
                      ## check whether all the linkage groups have been removed
                      if(length(private$LG) == 0){
                        stop("Linakge groups have not been formed. Use the '$createLG' function to create linkage groups.")
                        private$LG <- NULL
                      }
                      return(invisible(NULL))
                    },  
                    ## function for merging linkage groups
                    mergeLG = function(LG){
                      ## checks
                      if(is.null(private$LG)) 
                        stop("Linakge groups have not been formed. Use the '$createLG' function to create linkage groups.")
                      nLG <- length(private$LG)
                      ## check input
                      if(GUSbase::checkVector(LG, type="pos_integer", minv=1,maxv=nLG))
                        stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nLG))
                      LG <- unique(LG)
                      if(length(LG) == 1)
                        stop("Only one unique linkage group supplied. Need at least two linkage groups to merge")
                      ## merge LGs into the first
                      private$LG[[LG[1]]] <- unlist(private$LG[LG])
                      ## drop remaining groups
                      private$LG[LG[-1]] <- NULL
                      return(invisible(NULL))
                    },
                    ## function for masking SNPs
                    maskSNP = function(snps){
                      ## Inout checks
                      if( !is.vector(snps) || !is.numeric(snps) || any(is.na(snps)) || !all(snps == round(snps)) || any(snps < 1) || any(snps > private$nSnps) )
                        stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                      ## mask SNPs
                      private$masked[snps] <- TRUE
                      return(invisible(NULL))
                    },
                    ## function for unmasking SNPs
                    unmaskSNP = function(snps){
                      if( !is.vector(snps) || !is.numeric(snps) || any(is.na(snps)) || !all(snps == round(snps)) || any(snps < 1) || any(snps > private$nSnps) )
                        stop(paste0("Input must be a vector of indices between 1 and ", private$nSnps))
                      private$masked[snps] <- FALSE
                      return(invisible(NULL))
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
                                             private$config_orig[[1]], private$config_infer[[1]],
                                             private$group, private$group_infer,
                                             nClust, private$nInd, err=err)
                      else
                        stop("Multiple families have yet to be implemented")
                      #mat <- rf_2pt_multi(private$ref, private$alt,
                      #                    private$config,private$group, nClust, private$noFam)
                      ## Save the results to the object
                      private$rf <- mat$rf
                      private$LOD <- mat$LOD
                      return(invisible(NULL))
                    },
                    ## Function for creating linkage groups
                    createLG = function(LODthres=10, nComp=10, reset=FALSE){
                      ## Do some checks
                      if(is.null(private$rf) || is.null(private$LOD))
                        stop("Recombination fractions and LOD scores have not been computed.\nUse $rf_2pt() to compute the recombination fractions and LOD scores.")
                      if(!is.numeric(LODthres) || !is.vector(LODthres) || length(LODthres) != 1 || is.na(nComp) ||  LODthres < 0 || !is.finite(LODthres))
                        stop("The LOD threshold (argument 2) needs to be a finite numeric number.")
                      if(!is.numeric(nComp) || !is.vector(nComp) || length(nComp) != 1 || is.na(nComp) || nComp < 0 || !is.finite(nComp) ||
                         round(nComp) != nComp)
                        stop("The number of comparsion points (argument 3) needs to be a finite integer number.")
                      
                      ## reset
                      if(reset){
                        private$LG_map = NULL
                        private$para = NULL
                        private$summaryInfo = private$summaryInfo$data
                      }
                        
                      ## Create the groups
                      private$LG <- createLG(private$group, private$LOD, "both", LODthres, nComp, masked=private$masked)
                      ## display the results
                      self$print(what = "LG")
                      return(invisible(NULL))
                    },
                    ## Function for ordering linkage groups
                    orderLG = function(LG = NULL, mapfun = "haldane", weight="LOD2", ndim=30, spar=NULL, filename=NULL){
                      ## do some checks
                      if(is.null(private$LG))
                        stop("Linkage groups have not been formed. Use the '$createLG' function to create linkage groups.")
                      if(is.null(LG))
                        LG <- 1:length(private$LG)
                      else if(!is.vector(LG) && !is.numeric(LG) && any(LG < 0) && LG > length(private$LG))
                        stop("Invalid Linkage Group input")
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
                      nChr = length(LG)
                      ## order the linkage groups
                      for(lg in LG){
                        ind <- private$LG[[lg]]
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
                        graphics::image(private$rf[ind,ind][pcurve$ord,pcurve$ord], axes=F, col=grDevices::heat.colors(100))
                        if(is.null(filename)) graphics::par(temp_par)
                        else grDevices::dev.off()
                        ## Set the new order
                        private$LG[[lg]] <- ind[pcurve$ord]
                      }
                      ## Acknowledge paper we are using algorithm from
                      cat("Using the MDS scaling algorithm to order SNPs.\n",
                          "Please cite:\n",
                          "\tPreedy & Hackett (2016). Theor Appl Genet. 129:2117-2132\n",sep="")
                      return(invisible(NULL))
                    },
                    ## Function for plotting linkage groups
                    plotLG = function(LG=NULL, mat="rf", filename=NULL, interactive=TRUE, ...){
                      ## do some checks
                      if(!is.vector(mat) || !is.character(mat) || length(mat) != 1 || !(mat %in% c('rf','LOD')))
                        stop("Argument specifying which matrix to plot (argument 1) must be either 'rf' or 'LOD'")
                      if(!is.vector(interactive) || length(interactive) != 1 || !is.logical(interactive))
                        stop("Argument sepcifting whether to produce an interactive heatmap or not is invalid.")
                      if(!is.null(filename) && (!is.vector(filename) || !is.character(filename) || length(filename) != 1))
                        stop("Specified filename is invalid")
                      plotly.arg <- list(...)
                      if(is.null(plotly.arg$selfcontained))
                        plotly.arg$selfcontained = FALSE
                      
                      if(private$noFam == 1){
                        ## Work out which LGs list to use
                        if(is.null(private$LG))
                          stop("No linkage groups are available to be plotted. Please use the $createLG function to create some linkage groups")
                        else
                          LGlist <- c(private$LG)
                        ## Work out if we want a subset of the LGs
                        if(!is.null(LG)){
                          if(GUSbase::checkVector(LG, type="pos_integer", minv=1, maxv=length(LGlist)))
                            stop(paste0("At least one linkage group number does not exist. Indices must be between 1 and",length(LGlist),"\n"))
                          LGlist <- LGlist[LG]
                          names(LGlist) <- LG
                        }
                        else
                          names(LGlist) <- 1:length(LGlist)
                        
                        ## Sort out the matrix
                        if(mat == "rf"){
                          temprf <- private$rf
                          heatmapcol = grDevices::heat.colors(100)
                          zrange = c(0,0.5)
                        }
                        else if(mat == "LOD"){
                          temprf <- private$LOD
                          heatmapcol = rev(grDevices::heat.colors(100))
                          zrange = c(0,min(50,max(private$LOD)))
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
                          # suppress warnings  
                          storeWarn <- getOption("warn")
                          options(warn = -1) 
                          ## produce the plotly plot
                          if(length(which(b_indx)) == 0){
                            p <- plotly::plot_ly(z=temprf, type="heatmap", showscale=F, hoverinfo="text",
                                                 text=hovertext, colors=heatmapcol) %>%
                              plotly::layout(margin=list(l=0,r=0,t=0,b=0), xaxis=ax, yaxis=ax)
                          }
                          else{
                            p <- plotly::plot_ly(z=temprf, type="heatmap", showscale=F, hoverinfo="text",
                                                 text=hovertext, colors=heatmapcol) %>% 
                              plotly::add_segments(x=which(b_indx)-1,xend=which(b_indx)-1,y=0,yend=nn, line=list(color="black"),  showlegend=F) %>%
                              plotly::add_segments(y=which(b_indx)-1,yend=which(b_indx)-1,x=0,xend=nn, line=list(color="black"),  showlegend=F) %>%
                              plotly::layout(margin=list(l=0,r=0,t=0,b=0), xaxis=ax, yaxis=ax)
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
                            temp_par = graphics::par(no.readonly = TRUE)
                          if(mat == "LOD") {
                            highrf = which(temprf > max(zrange))
                            if(length(highrf) > 0) temprf[highrf] = max(zrange)
                          }
                          graphics::par(mar=rep(0,4), oma=c(0,0,0,0), mfrow=c(1,1), xaxt='n', yaxt='n', bty='n',ann=F)
                          graphics::image(temprf, x=1:nn, y=1:nn, zlim=zrange, col=heatmapcol)
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
                    plotChr = function(chrom = NULL, mat=c("rf"), filename=NULL, chrS=2, lmai=2){
                      ## do some checks
                      if(!is.vector(mat) || !is.character(mat) || length(mat) != 1 || !(mat %in% c('rf','LOD')))
                        stop("Argument specifying which matrix to plot (argument 1) must be either 'rf' or 'LOD'")
                      nChrom = length(unique(private$chrom))
                      if(is.null(chrom)){
                        chrom <- 1:nChrom
                      } else if(GUSbase::checkVector(chrom, type = "pos_integer", maxv=nChrom))
                        stop(paste0("Invalid chromosome number (first argument). Must be an integer number between 0 and ",nChrom,".\n"))
                      
                      if(is.null(filename))
                        temp_par <- par(no.readonly = TRUE) # save the current plot margins
                      
                      ## workout which SNPs to plot
                      if(private$noFam == 1){
                        ## workout the indices for the chromosomes
                        names = unique(private$chrom)
                        LGlist <- sapply(names, function(x) which((private$chrom == x) & !private$masked), simplify=F)
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
                      if(is.null(private$para))
                        stop("There are no linkage maps availabe to plot. Use the '$computeMap' function to compute linkage maps.")
                      nLGs = length(private$para[[1]])
                      if(is.null(LG))
                        LG = 1:nLGs
                      else if(GUSbase::checkVector(LG, type="pos_integer", minv=1, maxv=nLGs))
                        stop(paste0("The LG indices must be an integer vector between 1 and the number of linkage groups which is ",nLGs))
                      LG = LG[which(!unlist(lapply(private$para$rf[LG], function(x) any(is.na(x)))))]
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
                          col = rgb(t(col2rgb(col)), max=255)
                        if(length(col) == 1)
                          col <- rep(col,nLGs)
                        else if(length(col) != nLGs)
                          stop("Number of colors specified does not equal the number of linkage groups")
                      }
                      temp_par = par(no.readonly = TRUE) # save the current plot margins
                      
                      ellipseEq_pos = function(x) c2 + sqrt(b^2*(1-round((x-c1)^2/a^2,7)))
                      ellipseEq_neg = function(x) c2 - sqrt(b^2*(1-round((x-c1)^2/a^2,7)))
                      mapDist = lapply(private$para$rf[LG],mfun, fun=fun, centi=TRUE)
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
                        stop("There are no combined linkage groups availabe. Use the '$createLG' function to create linkage groups.")
                      ## work out the LG and original orders
                      LGorder <- unlist(private$LG)
                      orgOrder <- sort(LGorder)
                      ## determine the breask on the plot
                      temp <- cumsum(unlist(lapply(private$LG, length)))
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
                      graphics::mtext(text = 1:length(private$LG), side = 1, 
                                      at = apply(cbind(c(min(LGorder),LGbreaks),c(LGbreaks,max(LGorder))),1,mean))
                      graphics::par(temp_par) # reset the plot margins
                    }, 
                    ## Function for computing the rf's for each chromosome 
                    computeMap = function(chrom=NULL, init_r=0.001, ep=0.001, method="EM", err=TRUE, multiErr=FALSE, 
                                          mapped=TRUE, nThreads=1, inferOPGP=TRUE, rfthres=0.1, ...){
                      ## do some checks
                      if( !is.null(init_r) & !is.numeric(init_r) )
                        stop("Starting values for the recombination fraction needs to be a numeric vector or integer or a NULL object")
                      if( (length(ep) != 1 || !is.numeric(ep) || (ep <= 0 | ep >= 1)) )
                        stop("Value for the error parameters needs to be a single numeric value in the interval (0,1) or a NULL object")
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
                              tempOPGP <- c(tempOPGP,list(as.integer(infer_OPGP_FS(ref_temp[[fam]],alt_temp[[fam]],private$config[[fam]][indx_chrom], method="EM", nThreads=1, ...))))
                            }
                            private$para$OPGP[i] <- tempOPGP
                          }
                          ## estimate the rf's
                          MLE <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                           sexSpec=FALSE, seqErr=err, method=method, nThreads=nThreads, multiErr=multiErr)
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
                        if(length(private$LG) == 0)
                          stop("No linkage groups have been formed.")
                        nchrom <- length(private$LG)
                        if(is.null(chrom)) # compute distances for all chromosomes
                          chrom <- 1:nchrom
                        else if(GUSbase::checkVector(chrom, type="pos_integer", minv=0, maxv=nchrom))
                          stop("Linkage group input must be an integer vector between 1 and the number of linkage groups")
                        ## Compute rf's
                        cat("Computing recombination fractions:\n")
                        if(is.null(private$para) || length(private$para$rf) != length(private$LG))
                          private$para <- list(OPGP=replicate(nchrom, NA, simplify=FALSE),
                                               rf=replicate(nchrom, NA, simplify=FALSE),
                                               ep=replicate(nchrom, NA, simplify=FALSE),
                                               loglik=replicate(nchrom, NA, simplify=FALSE))
                        if(is.null(private$LG_map) || length(private$LG_map) != nchrom)
                          private$LG_map <- vector(mode="list", length = nchrom)
                        for(i in chrom){
                          cat("Linkage Group: ",i,"\n")
                          indx_chrom <- private$LG[[i]]
                          ref_temp <- lapply(private$ref, function(x) x[,indx_chrom])
                          alt_temp <- lapply(private$alt, function(x) x[,indx_chrom])
                          ## estimate OPGP's
                          curOPGP <- private$para$OPGP[i]
                          if(inferOPGP || !is.list(curOPGP) || length(curOPGP[[1]]) != length(indx_chrom)){
                            tempOPGP <- list()
                            for(fam in 1:private$noFam){
                              tempOPGP <- c(tempOPGP,list(as.integer(infer_OPGP_FS(ref_temp[[fam]],alt_temp[[fam]],private$config[[fam]][indx_chrom], method="EM", nThreads=1, ...))))
                            }
                            private$para$OPGP[i] <- tempOPGP
                          }
                          ## estimate the rf's
                          if(!is.null(method))
                            MLE <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                             sexSpec=FALSE, seqErr=err, method=method, nThreads=nThreads, multiErr=multiErr, ...)
                          else{
                            MLE_init <- rf_est_FS(init_r=init_r, ep=ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                                  sexSpec=FALSE, seqErr=err, method="EM", nThreads=nThreads, multiErr=multiErr, maxit=20)
                            MLE <- rf_est_FS(init_r=MLE_init$rf, ep=MLE_init$ep, ref=ref_temp, alt=alt_temp, OPGP=private$para$OPGP[i],
                                             sexSpec=FALSE, seqErr=err, method="optim", nThreads=nThreads, multiErr=multiErr, ...)
                          }
                          private$para$rf[i]   <- list(MLE$rf)
                          private$para$ep[i]       <- list(MLE$ep)
                          private$para$loglik[i]   <- list(MLE$loglik)
                          private$LG_map[[i]]      <- private$LG[[i]]
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
                        ## Create summary of results:
                        BI <- unlist(lapply(private$LG_map,length))
                        tab <- cbind(LG=1:nchrom,BI)
                        DIST <- extendVec(unlist(lapply(private$para$rf, function(x) round(sum(mfun(x, fun="haldane", centiM = T)),2))), nrow(tab))
                        ERR <- extendVec(round(unlist(lapply(private$para$ep,mean)),5), nrow(tab))
                        tab <- cbind(tab,DIST,ERR)
                        tab <- stats::addmargins(tab, margin = 1, FUN=function(x) sum(x, na.rm=T))
                        rownames(tab) <- NULL
                        tab[nrow(tab), 1] <- "TOTAL"
                        tab[nrow(tab), ncol(tab)] = ""
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
                          else if(!is.null(private$LG)) what = "LG"
                        }
                        ## Find out which LGs to write to file
                        if(what == "map"){
                          if(!is.vector(inferGeno) || length(inferGeno) != 1 || !is.logical(inferGeno) || is.na(inferGeno))
                            stop("Argument `inferGeno` is invalid. Should be a logical value")
                          if(is.null(private$LG_map) && is.null(private$para))
                            stop("No linkage maps available to write to file.")
                          if(is.null(LG)){
                            LGlist <- private$LG_map
                            LG <- 1:length(private$LG_map)
                          } else if(GUSbase::checkVector(LG, type="pos_integer", minv=1, maxv=length(private$LG_map)))
                            stop("LGs to write to file is invalid.")
                          else
                            LGlist <- private$LG_map[LG]
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
                          segType <- rep("BI",length(unlist(LGlist)))
                          temp <- private$ref[[1]][,unlist(LGlist)] > 0 | private$alt[[1]][,unlist(LGlist)] > 0
                          callrate <- colMeans(temp)
                          rf <- unlist(lapply(private$para$rf[LG], function(y) {if(!is.na(y[1])) return(format(round(cumsum(c(0,y)),6),digits=6,scientific=F))} ))
                          #rf_m <- unlist(lapply(private$para$rf_m[LG], function(y) {if(!is.na(y[1])) return(format(round(cumsum(c(0,y)),6),digits=6,scientific=F))} ))
                          ep <- format(round(unlist(sapply(1:length(LG), function(x) rep(private$para$ep[[x]], length.out=length(LGlist[[x]])))),8),digits=8,scientific=F)
                          ## add geno information if required
                          if(inferGeno){
                            temp = do.call("rbind", geno)
                            temp2 = split(temp, rep(1:ncol(temp), each = nrow(temp)))
                            names(temp2) = paste0(c("MAT_","PAT_"), rep(private$indID[[1]], rep(2,nInd)))
                            out <- c(list(LG=LGno, LG_POS=LGindx, CHROM=chrom, POS=pos, TYPE=segType, RF=rf,
                                          ERR=ep, MEAN_DEPTH=depth, CALLRATE=callrate), temp2)
                          } else 
                            out <- list(LG=LGno, LG_POS=LGindx, CHROM=chrom, POS=pos, TYPE=segType, RF=rf,
                                        ERR=ep, MEAN_DEPTH=depth, CALLRATE=callrate)
                          ## write out file
                          data.table::fwrite(out, file = outfile, sep="\t", nThread = 1)  
                          cat(paste0("Linkage analysis results written to file:\nFilename:\t",filename)) 
                          return(invisible(NULL))
                        } else if(what == "LG"){
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
                          segType <- rep("BI",length(unlist(LGlist)))
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
                    config       = NULL,
                    group        = NULL,
                    masked       = NULL,
                    noFam        = NULL,
                    rf           = NULL,
                    LOD          = NULL,
                    famInfo      = NULL,
                    para         = NULL,
                    LG           = NULL,
                    LG_map       = NULL,
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