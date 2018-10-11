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
#' FS method: Heatmap of 2-point rf and LOD estimates from linkage group order
#' 
#' Method for plotting 2-point recombination fraction and LOD score estimates when the SNPs are ordered according to the linkage groups.
#'
#' This function plots the heatmap of the matrix of 2-point recombination fraction estimates (or 2-point LOD scores)
#' as computed from the \code{\link{$rf_2pt}} function. 
#' 
#' The interactive plot (\code{interactive = TRUE}) is produced via the \code{\link[plotly]{plot_ly}} function. Note, however,
#' that the interactive plot only works for a relatively small number of SNPs (about 1000 is maximum).
#' Nevertheless, one can still plot each linkage group (using the \code{LG} argument) provided there are 
#' not too many SNPs in the linkage group.
#' 
#' When \code{what = "LG"}, the pseudo-testcross linkage groups created via the \code{\link{$createLG}} function 
#' will be plotted. On the other hand, if \code{what = "LG_BI"}, then the combined linkage groups created from 
#' the \code{\link{$addBIsnps}} function will be plotted. When \code{what = NULL}, the combined linkage groups will
#' be plotted if available, otherwise the pseudo-testcross linkage groups will be plotted.
#' 
#' @usage
#' FSobj$plotLG(parent  = "maternal", LG=NULL, mat="rf", filename=NULL, interactive=TRUE, what = NULL)
#' 
#' @param parent Character value specifying whether the SNPs segreagting in the maternal parent should be 
#' plotted (\code{"maternal"}), or whether whether the SNPs segreagting in the paternal parent should be 
#' plotted (\code{"paternal"}), or whether all the SNPs should be plotted (\code{"both"}).
#' @param LG Integer vector giving the indices of the linkage groups to be plotted.
#' @param mat Charater value for the matrix to be plotted. \code{"rf"} plots the matrix of 2-point 
#' recombination fractions while \code{"LOD"} plots the matrix of 2-point LOD scores.
#' @param filename Character value giving the name of the file to save the plot to. If \code{NULL}, the plot is displayed 
#' in the graphics window and not saved to a file. 
#' @param interative Logical value. If \code{TRUE} then an interactive plot is produced, otherwise a standard
#' base R plot is used.
#' @param what Character vector specifying which list of linkage groups to plot. \code{"LG-pts"} is for 
#' the pseudo-testcross linkage groups and \code{"LG-comb"} is for the combined linkage groups.
#' 
#' @name $plotLG
#' @author Timothy P. Bilton
#' @seealso \code{\link{FS}}
#' @examples 
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(replicate(2, sample(c(1,2,4), size=30, replace=TRUE), simplify=FALSE))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50)
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt()
#' ## create paternal and maternal linkage groups
#' F1data$createLG()
#' 
#' ## Plot the linkage groups
#' F1data$plotLG()
#' 
#' ## Plot the linkage groups: suppress interactive plot 
#' F1data$plotLG(interactive = FALSE)
#' @aliases $plotLG

plotLG <- function(mat, LG, filename=NULL, names=NULL, chrS=2, lmai=2, chrom=T, type="rf"){
  
  if(length(LG) == 1){
    chrS = 0; lmai = 0
  }
  
  b <- ncol(mat) + 1
  if(chrom){
    chrom.ind <- unlist(lapply(LG, function(x) c(x,b)))[-length(unlist(LG))+length(LG)]
    chrom.ind <- chrom.ind[-length(chrom.ind)]
  }
  else
    chrom.ind <- 1:ncol(mat)
  
  ## Subset the matrix
  mat <- cbind(mat,rep(0.5,b-1))
  mat <- rbind(mat,rep(0.5,b))
  mat <- mat[chrom.ind,chrom.ind]
  ## work out where the breaks are
  breaks <- which(chrom.ind==b)
  npixels <- length(chrom.ind)
  if(chrom){
    if(!is.null(filename))
      grDevices::png(filename,width=npixels+72*lmai,height=npixels,res=72)
    graphics::par(mfrow=c(1,1), xaxt='n',yaxt='n',mai=c(0,lmai,0,0),bty='n',ann=F)
    ## what matrix to plot
    if(type == "rf")
      graphics::image(1:npixels,1:npixels,mat,zlim=c(0,0.5),col=grDevices::heat.colors(100))
    else if (type == "LOD")
      graphics::image(1:npixels,1:npixels,mat,zlim = c(0,50), 
            col=grDevices::colorRampPalette(rev(c("red","orange","yellow","white")), bias=5)(100))
    if(is.null(names))
      graphics::mtext(paste("LG",1:length(LG),"  "),
                  at=floor(apply(cbind(c(0,breaks),c(breaks,npixels)),1,stats::median)),side=2, line=0,cex=chrS,las=1)
    else
      graphics::mtext(names, at=floor(apply(cbind(c(0,breaks),c(breaks,npixels)),1,stats::median)),side=2, line=0,cex=chrS,las=1)
    graphics::abline(h=breaks)
    graphics::abline(v=breaks)
    if(!is.null(filename))
      grDevices::dev.off()
  }
  else{
    npixels <- length(chrom.ind)
    if(!is.null(filename))
      grDevices::png(filename,width=npixels,height=npixels)
    ## what matrix to plot
    if(type == "rf")
      graphics::image(1:npixels,1:npixels,mat,zlim=c(0,0.5),col=grDevices::heat.colors(100))
    else if (type == "LOD")
      graphics::image(1:npixels,1:npixels,mat,zlim = c(0,50), 
            col=grDevices::colorRampPalette(rev(c("red","orange","yellow","white")), bias=5)(100))
    graphics::abline(h=breaks)
    graphics::abline(v=breaks)
    if(!is.null(filename))
      grDevices::dev.off()
  }
}
