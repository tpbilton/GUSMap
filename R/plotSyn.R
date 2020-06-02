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
#' BC, FS and IC method: Produce synteny plot
#' 
#' Method for plotting the combined linkage group SNP ordering against the original SNP order
#' from the genomic assembly.
#' 
#' This function produces a scatter plot of the SNP order as obtained from the linkage analysis
#' against the SNP order given by the genomic assembly. 
#'
#' \code{BC} object: The SNP order from the linkage analysis is from the pseudo-testcross linkage groups with BI SNPs
#' constructed from the \code{\link{addBIsnps}} function and ordered via the \code{\link{$orderLG}} function.
#' 
#' \code{FS} object: The SNP order from the linkage analysis is from the combined linkage groups 
#' obtained from the \code{\link{$addBIsnps}} function and ordered via the \code{\link{$orderLG}} function.
#' 
#' \code{IC} object: The SNP order from the linkage analysis is from the linkage groups 
#' computed from the \code{\link{$computeLG}} function and ordered via the \code{\link{$orderLG}} function.
#' 
#' @usage
#' BCobj$plotSyn()
#' FSobj$plotSyn()
#' ICobj$plotSyn()
#' 
#' @name $plotSyn
#' @author Timothy P. Bilton
#' @seealso \code{\link{BC}}, \code{\link{FS}}, \code{\link{IC}}
#' @examples
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(replicate(2,sample(c(1,2,4), size=30, replace=TRUE), simplify=FALSE))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50)
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt()
#' ## create paternal and maternal linkage groups
#' F1data$createLG()
#' ## Add BI SNPs
#' F1data$addBIsnps()
#' ## order Markers
#' F1data$orderLG(ndim=5)
#' 
#' ## Plot synteny plot
#' F1data$plotSyn()
NULL
