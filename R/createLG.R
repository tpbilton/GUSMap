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
#' FS method: Create linkage groups
#' 
#' Method for creating maternal and paternal linakge groups using the computed LOD scores
#' 
#' Linkage groups are formed using the following algorithm
#' \enumerate{
#' \item The two unmapped SNPs (e.g. not in a linkage group) with the highest LOD score are grouped together to
#'  form first linkage group.
#' \item The average LOD score between each unmapped SNP and a specified number (\code{nComp}) of SNPs in the 
#' linkage group that have the highest LOD score with the unmapped SNP is computed. The unmapped SNP with the 
#' largest average LOD score is mapped to the linkage group if the average LOD score is above the LOD threshold
#' (given by \code{LODthres}).
#' \item Steps 1-3 are repeated until no more linkage groups can be formed or all the SNPs have been mapped 
#' to a linkage group.
#' }
#' 
#' The performance of the linakge group algorithm depends on the value of \code{LODthres} and \code{nComp}. The user is
#' advised to experiment with different values and examine the matrix of recombination fractions (using the \code{\link{$plotLG}}
#' function). 
#' 
#' Note: The LOD scores used in this function are computed using the {\code{\link{$rf_2pt}}} which needs to be
#' run beforehand.
#' 
#' @usage
#' FSobj$createLG(parent = "both", LODthres = 10, nComp = 10)
#' 
#' @param parent A character vector specifying whether to create linakge groups using the maternal-informative (MI)
#' SNPs \code{parent = "maternal"}, the paternal-informative (PI) SNPs \code{parent = "paternal"}, or MI and PI SNPs
#' \code{parent = "both"}.
#' @param LODthres A positive numeric value specifying the LOD threshold used to add SNPs to the linkage groups.
#' @param nComp A positive integer value specifying how many SNPs in the linakge group to compute the average LOD score 
#' with the unmapped SNP.
#' 
#' @name $createLG
#' @author Timothy P. Bilton
#' @seealso \code{\link{FS}}
#' @examples 
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(list(sample(c(1,2,4), size=30, replace=T)))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50)
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt(nClust=1)
#' 
#' ## create paternal and maternal linkage groups
#' F1data$createLG()
NULL

















createLG <- function(group, LOD, parent, LODthres, nComp, masked){

  ## Initialize the LG
  LG <- list()
  
  ## check that the parent argument is correct
  if(!is.character(parent) || length(parent) != 1 || 
     !(parent %in% c("maternal","paternal")))
    stop("parent argument is not a string of of length one or is incorrect:
         Please select one of the following:
         maternal: Only MI SNPs
         paternal: Only PI SNPs")
  
  if(parent == "maternal")
    unmapped <- sort(group$MI[which(group$MI %in% which(!masked))])
  else if(parent == "paternal")
    unmapped <- sort(group$PI[which(group$PI %in% which(!masked))])
  
  if(length(unmapped) < 2)
    stop("There are no SNPs available to create linkage groups with.")
  
  ## Run algorithm for generating the linkage groups
  finish = FALSE
  while(!finish){
    newLG <- unmapped[sort(which(LOD[unmapped,unmapped]==max(LOD[unmapped,unmapped],na.rm=T),arr.ind=T)[1,])]
    unmapped <- unmapped[-which(unmapped%in%newLG)]
    compLOD <- LOD[newLG,unmapped]
    meanLOD <- apply(compLOD,2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
    maxMeanLOD <- max(meanLOD)
    while(maxMeanLOD > LODthres){
      newSNP <- which(meanLOD == maxMeanLOD)
      newLG <- c(newLG,unmapped[newSNP])
      unmapped <- unmapped[-newSNP]
      if(length(unmapped) > 0){
        compLOD <- matrix(LOD[newLG,unmapped],ncol=length(unmapped))
        meanLOD <- apply(compLOD,2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
        maxMeanLOD <- max(meanLOD)
      }
      else {
        finish = TRUE
        maxMeanLOD = -999
        next
      }
    }
    meanLOD <- apply(matrix(LOD[unmapped,unmapped], nrow=length(unmapped), ncol=length(unmapped)),
                     2,function(x) mean(sort(x,decreasing = T)[1:nComp], na.rm=T))
    finish <- !any(meanLOD > LODthres)
    LG <- c(LG,list(newLG))
  }
  return(LG)
}
