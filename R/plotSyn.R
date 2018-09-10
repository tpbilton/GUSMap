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
#' FS method: Produce synteny plot
#' 
#' Method for plotting the combined linkage group SNP ordering against the original SNP order
#' from the genomic assembly.
#' 
#' This function produces a scatter plot of the SNP order as obtained from the linkage analysis
#' against the SNP order given by the genomic assembly. Note that the SNP order from the linkage analysis
#' is from the combined linkage groups obtained using \code{\link{$plotSyn}} and ordered via \code{\link{$orderLG}}.
#' 
#' @usage
#' FSobj$plotSyn()
#' 
#' @name $plotSyn
#' @author Timothy P. Bilton
#' @seealso \code{\link{FS}}
#' @examples
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(replicate(2,sample(c(1,2,4), size=30, replace=T), simplify=FALSE))
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
