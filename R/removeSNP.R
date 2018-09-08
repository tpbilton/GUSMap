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
#' FS method: Remove SNP from Linkage Group
#' 
#' Method for removing SNPs from linkage groups in an FS object.
#' 
#' For a linkage analysis in GUSMap, there may be a need to delete linkage groups from the
#' list of created linkage groups. The indices of the linkage groups corresponds to the number
#' given in the output for the FS object. More than one linkage group can be deleted at once.
#' 
#' Note: When the linkage group(s) are deleted, the indices of the remaining linkage groups 
#' changes. Print the FS object to investiage what these after running this function.
#' 
#' @usage
#' FSobj$removeSNP(snps)
#' 
#' @param snps An integer vector giving the indices of the SNP to be removed.
#' @author Timothy P. Bilton
#' @name $removeSNP
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
#' ## Remove SNPs 3,7,15,20 from the linkage groups
#' F1data$removeSNP(snps = c(3,7,15,20))
NULL