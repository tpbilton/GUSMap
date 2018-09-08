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
#' FS method: Write linkage groups/maps to file
#'
#' Method for writing the results from the linkage analysis in GUSMap to a file.
#' 
#' When \code{what = "LG"}, results for the linkage groups from the paternal and maternal linkage 
#' group list created via the \code{\link{$createLG}} function will be return. On the other hand, if 
#' \code{where = "LG_BI"}, linkage groups from the combined linkage group list 
#' created from the \code{\link{$addBIsnps}} function (e.g., the linkage groups which include the BI SNPs) will be returned
#' and the linkage maps (if they have been produced using the \code{\link{$computeMap}} function) will also be returned.
#' When \code{where = NULL}, the results from the combined linkage groups will be returned if available,
#' otherwise the result for the maternal and paternal linkage group list wil be returned.
#' 
#' @usage
#' FSobj$writeLM(file, direct = "./", LG = NULL)
#' 
#' @name $writeLM
#' @param file Character value giving the name of the file to write to.
#' @param direct Character for the directory to write the file to (relative to the current working working directory)
#' @param LG Integer vector giving the indices of the linkage groups to write. If \code{NULL}, the results for all the 
#' linkage groups are returned.
#' @return 
#' The function returns a csv file with the following columns:
#' \describe{
#' \item{LG}{The linkage number as given in the output from the FS object.}
#' \item{LGPOS}{The position of the SNPs across the linkage group. Starts at 1 and
#' increases by one for each SNP.}
#' \item{CHROM}{The chromosome name the SNP was located on in the original genomic assembly.}
#' \item{POS}{The position (in base pairs) the SNP was located on in the original genomic assembly.}
#' \item{TYPE}{The segregation type of SNP. Either BI (both-informative), 
#' MI (maternal-informative) or PI (paternal-informative).}
#' \item{RF_PAT}{The sum of the recombination fraction estimates for the paternal meioses.}
#' \item{RF_MAT}{The sum of the recombination fraction estimates for the maternal meioses.}
#' \item{ERR}{THe sequencing error estimate for the SNP.}
#' \item{MEAN_DEPTH}{The mean read depth of the SNP.}
#' \item{CALLRATE}{The proportion of individuals in which there is at least one read across the SNP.}
#' }
#' If no linkage maps have been computed (using the \code{\link{$computeMap}}), then the
#' columns RF_PAT, RF_MAT and ERR are not returned.
#' 
#' @author Timothy P. Bilton
#' @seealso \code{\link{FS}}
#' @examples
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(list(sample(c(1,2,4), size=30, replace=T)))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50)
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt(nClust=1)
#' ## create and order linkage groups
#' F1data$createLG()
#' F1data$addBIsnps()
#' F1data$orderLG(ndim=5)
#' F1data$computeMap()
#' 
#' ## return the results to a file
#' F1data$writeLM(file = "test")
NULL