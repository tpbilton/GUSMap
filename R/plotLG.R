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
#' FS method: Plot linkage groups
#' 
#' Plot linkage groups
#' 
#' 
#' @name $plotLG
#' @author Timothy P. Bilton
#' @seealso \code{\link{FS}}
#' @examples 
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(list(sample(c(1,2,4), size=10, replace=T)), list(sample(c(1,2,4), size=10, replace=T)))
#' F1data <- simFS(0.01, config=config, meanDepth=10, 
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
#' @aliases NULL

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
      png(filename,width=npixels+72*lmai,height=npixels,res=72)
    par(xaxt='n',yaxt='n',mai=c(0,lmai,0,0),bty='n',ann=F)
    ## what matrix to plot
    if(type == "rf")
      image(1:npixels,1:npixels,mat,zlim=c(0,0.5),col=heat.colors(100))
    else if (type == "LOD")
      image(1:npixels,1:npixels,mat,zlim = c(0,50), 
            col=colorRampPalette(rev(c("red","orange","yellow","white")), bias=5)(100))
    if(is.null(names))
      mtext(paste("LG",1:length(LG),"  "),
                  at=floor(apply(cbind(c(0,breaks),c(breaks,npixels)),1,median)),side=2, line=0,cex=chrS,las=1)
    else
      mtext(names, at=floor(apply(cbind(c(0,breaks),c(breaks,npixels)),1,median)),side=2, line=0,cex=chrS,las=1)
    abline(h=breaks)
    abline(v=breaks)
    if(!is.null(filename))
      dev.off()
  }
  else{
    npixels <- length(chrom.ind)
    if(!is.null(filename))
      png(filename,width=npixels,height=npixels)
    ## what matrix to plot
    if(type == "rf")
      image(1:npixels,1:npixels,mat,zlim=c(0,0.5),col=heat.colors(100))
    else if (type == "LOD")
      image(1:npixels,1:npixels,mat,zlim = c(0,50), 
            col=colorRampPalette(rev(c("red","orange","yellow","white")), bias=5)(100))
    abline(h=breaks)
    abline(v=breaks)
    if(!is.null(filename))
      dev.off()
  }
}
