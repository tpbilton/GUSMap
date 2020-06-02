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
#' BC and FS method: Add MI, PI and SI SNPs to existing linkage groups
#' 
#' Method for adding maternal-informative (MI), paternal-informative (PI) and semi-informative (SI) SNPs 
#' to the paternal and maternal (pseudo-testcross) linkage groups.
#' 
#' Each MI, PI or SI SNP is mapped to the maternal and paternal linkage groups using the following
#' algorithm:
#' \enumerate{
#' \item For each pseudo-testcross linkage group, compute the average LOD score between the unmapped SNP and a 
#' specified number (\code{nComp}) of SNPs in the pseudo-testcross linkage group that have the highest LOD score 
#' with the unmapped SNP.
#' \item Map the unmapped SNP to the pseudo-testcross linkage group with the highest average LOD score if:
#' \itemize{
#'   \item The largest average LOD score is greater than the LOD threshold (\code{LODthres}).
#'   \item All the other average LOD scores associated with the other linkage groups are less than the LOD threshold.
#'   }
#' }
#' 
#' @usage
#' BCobj$addSNPs(LODthres = 10, nComp = 30)
#' FSobj$addSNPs(LODthres = 10, nComp = 30)
#' 
#' @param LODthres A positive numeric value specifying the LOD threshold used to add SNPs to the linkage groups.
#' @param nComp A positive integer value specifying how many SNPs in the linkage group to compute the average LOD score 
#' with the unmapped SNP.
#' 
#' @author Timothy P. Bilton
#' @name $addSNPs
#' @seealso \code{\link{BC}}, \code{\link{FS}}
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
#' ## Add more SNPs to the linkage groups
#' F1data$addSNPs()
NULL


