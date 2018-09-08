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
#'
#' FS method: Remove Linkage Groups
#' 
#' Method for removing linkage group(s) from the list of linkage groups contained
#' in an FS object.
#' 
#' In a linkage analysis in GUSMap, there may be need to remove linkage groups from the 
#' list of created linkage groups. The indices of the linkage groups corresponds to the number
#' given in the output for the FS object. 
#' 
#' When \code{where = "LG"}, linkage groups will be removed from the maternal and paternal linkage 
#' group list created via the \code{\link{$createLG}} function. On the other hand, if 
#' \code{where = "LG_BI"}, then linkage groups will be removed from the combined linkage group list 
#' created from the \code{\link{$addBIsnps}} function (e.g., the linkage group list which include the BI SNPs).
#' When \code{where = NULL}, linkage groups will be removed from the combined linkage groups if available,
#' otherwise they will be removed from the maternal and paternal linkage group list.
#' 
#' @usage
#' FSobj$removeLG(LG, what = NULL)
#' 
#' @param LG An integer vector specifying the number of the linkage groups to be removed.
#' @param where Character vector specifying which list of linkage groups to remove linkage groups from. \code{"LG"} is for 
#' the maternal and paternal linkage groups and \code{"LG_BI"} is for the combined linkage group list (see details).
#' @name $removeLG
#' @author Timothy P. Bilton
#' @seealso \code{\link{FS}}
#' @examples 
#' ## simulate sequencing data
#' config <- list(replicate(2,sample(c(1,2,4), size=30, replace=T), simplify = FALSE))
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