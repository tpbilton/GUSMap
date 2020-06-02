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
#' BC, FS and IC method: Remove linkage groups
#' 
#' Method for removing linkage group(s) from the list of linkage groups contained
#' in a BC, FS or IC object.
#' 
#' For a linkage analysis in GUSMap, there may be the need to remove linkage groups from the 
#' list of created linkage groups. The indices of the linkage groups corresponds to the number
#' given in the output for the BC, FS or IC object.
#'
#' \code{BC} object: When \code{where = "LG-pts"}, linkage groups will be removed from the set of pseudo-testcross 
#' linkage groups created via the \code{\link{$createLG}} function. On the other hand, if 
#' \code{where = "LG-bi"}, then the pseudo-testcross linkage groups with BI SNPs added created
#' from the \code{\link{$addBIsnps}} function will be removed.
#' When \code{where = NULL}, linkage groups will be removed from the set of pseudo-testcross linkage groups with BI SNPs added if available,
#' otherwise they will be removed from the set of pseudo-testcross linkage groups.
#' 
#' \code{FS} object: When \code{where = "LG-pts"}, linkage groups will be removed from the set of pseudo-testcross 
#' linkage groups created via the \code{\link{$createLG}} function. On the other hand, if 
#' \code{where = "LG-comb"}, then linkage groups from the set of combined linkage groups 
#' created from the \code{\link{$addBIsnps}} function will be removed.
#' When \code{where = NULL}, linkage groups will be removed from the set of combined linkage groups if available,
#' otherwise they will be removed from the set of pseudo-testcross linkage groups.
#' 
#' @usage
#' BCobj$removeLG(LG, where = NULL)
#' FSobj$removeLG(LG, where = NULL)
#' ICobj$removeLG(LG)
#' 
#' @param LG An integer vector specifying the number of the linkage groups to be removed.
#' @param where Character vector specifying which list of linkage groups to remove linkage groups from. \code{"LG-pts"} is for 
#' the maternal and paternal pseudo-testcross linkage groups, \code{"LG-bi"} is for the pseudo-testcross linkage groups with BI SNPs
#' added (\code{BC} objects only) and \code{"LG-comb"} is for the combined linkage group list (\code{FS} objects only).
#' @name $removeLG
#' @author Timothy P. Bilton
#' @seealso \code{\link{BC}}, \code{\link{FS}}, \code{\link{IC}}
#' @examples 
#' ## simulate sequencing data
#' config <- list(replicate(2,sample(c(1,2,4), size=30, replace=TRUE), simplify = FALSE))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50)
#' 
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt(nClust=1)
#' ## create paternal and maternal linkage groups
#' F1data$createLG()
#' 
#' ## remove linkage groups 2 and 4
#' F1data$removeLG(LG = c(2,4))
NULL