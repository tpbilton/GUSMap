##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2020 Timothy P. Bilton <timothy.biltn@agresearch.co.nz>
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
#' FS method: Add BI SNPs to existing linkage groups
#' 
#' Method for adding both-informative (BI) SNPs to paternal and maternal linkage groups
#' and to combine maternal and paternal linkage groups.
#' 
#' Each both-informative (BI) SNP is mapped to the maternal and paternal linkage groups using the following
#' algorithm:
#' \enumerate{
#' \item For each MI (or PI) linkage group, compute the average LOD score between the unmapped BI SNP and a 
#' specified number (\code{nComp}) of SNPs in the MI (or PI) linkage group that have the highest LOD score 
#' with the unmapped BI SNP.
#' \item Map the unmapped BI SNP to linkage group with the highest average LOD score if:
#' \itemize{
#'   \item The largest average LOD score is greater than the LOD threshold (\code{LODthres}).
#'   \item All the other average LOD scores are less than the LOD threshold.
#'   }
#' }
#' Each BI SNP is mapped to a MI linkage group and then independently to a PI linkage group using the above algorithm. The 
#' MI and PI linkage groups are then merged together using the following algorithm:
#' \enumerate{
#' \item Compute the number of BI SNPs that mapped to the same pair of MI and PI linkage groups.
#' \item Merge the pair of MI and PI linkage groups that has the largest number of BI SNPs mapped to them,
#' provided that both linkage groups have not already been merged.
#' \item Contiune step 2 until either every MI linkage group has been merged with a PI linkage group, or 
#' every PI linkage group has been merged with a MI linkage group, or until no more MI or PI linkage groups can 
#' can be merged to another non-merged linkage group.
#' \item Merge the remaining (unmerged) MI or PI linkage groups to the combined linkage group with the most BI SNPs 
#' in common (provided there are any).
#' \item Contiune step 4 until all linkage groups have been merged or until no more linkage groups can be merged.
#' }
#' The performance of the linkage group algorithm depends on the value of \code{LODthres} and \code{nComp}. The user is
#' advised to experiment with different values and examine the matrix of recombination fractions (using the \code{\link{$plotLG}}
#' function). 
#' 
#' Note: The linkage groups produced from this function are referred to as the "combined linkage groups" and are stored separately
#' to the "pseudo-testcross linkage groups" produced by \code{\link{$createLG}} function. This means that combined linkage groups
#' can produced at any time from the pseudo-testcross linkage groups, even after major edits to the combined linkage groups
#' have been made.
#' 
#' @usage
#' FSobj$addBIsnps(LODthres = 10, nComp = 30)
#' 
#' @param LODthres A positive numeric value specifying the LOD threshold used to add SNPs to the linkage groups.
#' @param nComp A positive integer value specifying how many SNPs in the linkage group to compute the average LOD score 
#' with the unmapped SNP.
#' 
#' @author Timothy P. Bilton
#' @name $addBIsnps
#' @seealso \code{\link{FS}}
#' @examples 
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(list(sample(c(1,2,4), size=30, replace=TRUE)))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd = 50)
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt()
#' ## create paternal and maternal linkage groups
#' F1data$createLG()
#' 
#' ## Add the BI SNPs
#' F1data$addBIsnps()
NULL