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
#' FS method: Mask and unmask SNPs 
#' 
#' Method for masking and unmasking SNPs in an FS object. 
#' 
#' In the linkage mapping process, some SNPs may some themselves to be probmatic and there 
#' may be the need exclude these SNPs in forming linkage groups etc. This can be achieved using
#' the \code{\link{$maskSNP}} function, which masks SNPs. 
#' Once a SNP is masked, it will be ignored in the functions \code{\link{$createLG}},
#' \code{\link{$addBIsnps}}, \code{\link{$computeMap}}. Note, however, that if a SNP is
#' already included in a linkage group, masking a SNP will not remove it (To remove SNPs from
#' a linkage group, one needs to use the \code{\link{$removeSNP}} function). 
#' 
#' The \code{\link{$unmaskSNP}} function just unmasks any SNPs which have been masked using 
#' \code{\link{$maskSNP}}.
#' 
#' @usage
#' FS$obj$maskSNP(snps)
#' 
#' FS$obj$unmaskSNP(snps)
#' 
#' @param snps An integer vector giving the index number of the SNP in the dataset.
#' 
#' @name $maskSNP
#' @aliases $unmaskSNP
#' @author Timothy P. Bilton
#' @seealso \code{\link{FS}}
#' @examples
#' ## simulate some sequencing data
#' set.seed(6745)
#' config <- list(list(sample(c(1,2,4), size=30, replace=T)))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50)
#'
#' ## Mask SNPs 1,5,8,11,13
#' F1data$maskSNP(c(1,5,8,11))
#' 
#' ## Unmask SNPs 5 and 8
#' F1data$unmaskSNP(c(5,8))
#'
NULL