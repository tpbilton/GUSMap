##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2020 Timothy P. Bilton <timothy.bilton@agresearch.co.nz>
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
#'
#' BC and FS method: Remove SNP from Linkage Group
#' 
#' Method for removing SNPs from linkage groups in an BC or FS object.
#' 
#' For a linkage analysis in GUSMap, there may be a need to delete linkage groups from the
#' list of created linkage groups. The indices of the linkage groups corresponds to the number
#' given in the output for the FS object. More than one linkage group can be deleted at once.
#' 
#' Note: When the linkage group(s) are deleted, the indices of the remaining linkage groups 
#' changes. Print the BC or FS object to investiage what these are after running this function.
#'
#' \code{BC} obsject: When \code{where = "LG-pts"}, SNPs will be removed from the set of pseudo-testcross 
#' linkage groups created via the \code{\link{$createLG}} function. On the other hand, if 
#' \code{where = "LG-bi"}, then SNPs from the set of pseudo-testcross linkage groups with BI SNPs added
#' created from the \code{\link{$addBIsnps}} function will be removed.
#' When \code{where = NULL}, SNPs will be removed from the set of pseudo-testcross linkage groups with BI SNPs added if available,
#' otherwise they will be removed from the set of pseudo-testcross linkage groups.
#'
#' \code{FS} obsject: When \code{where = "LG-pts"}, SNPs will be removed from the set of pseudo-testcross 
#' linkage groups created via the \code{\link{$createLG}} function. On the other hand, if 
#' \code{where = "LG-comb"}, then SNPs from the set of combined linkage groups 
#' created from the \code{\link{$addBIsnps}} function will be removed.
#' When \code{where = NULL}, SNPs will be removed from the set of combined linkage groups if available,
#' otherwise they will be removed from the set of pseudo-testcross linkage groups.
#' 
#' @usage
#' BCobj$removeSNP(snps, where = NULL)
#' FSobj$removeSNP(snps, where = NULL)
#' 
#' @param snps An integer vector giving the indices of the SNP to be removed.
#' @param where Character vector specifying which set of linkage groups to remove SNPs from (see details).
#' @author Timothy P. Bilton
#' @name $removeSNP
#' @seealso \code{\link{BC}}, \code{\link{FS}}
#' @examples 
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(list(sample(c(1,2,4), size=20, replace=TRUE)))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50) 
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt()
#' ## create paternal and maternal linkage groups
#' F1data$createLG()
#' 
#' ## Remove SNPs 3,7,15,20 from the linkage groups
#' F1data$removeSNP(snps = c(3,7,15,20))
NULL