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
#' BC, FS and IC method: Write linkage groups/maps to file
#'
#' Method for writing results from a linkage analysis in GUSMap to a file.
#' 
#' Linkage groups and linkage maps (if computed) are written to a file. 
#' When \code{what = NULL}, the linkage mapping results computed using the \code{\link{$computeMap}} function are returned
#' if they are available, otherwise the linkage groups with BI SNPs information is returned. 
#' 
#' If \code{inferGeno=TRUE}, then the maternal and paternal haplotypes of the progeny in the population are returned. These haplotypes
#' are computed using the Viterbi algorithm with the parameters from the final map.  
#' 
#' @usage
#' BCobj$writeLM(file, direct = "./", LG = NULL, what = NULL, inferGeno = TRUE)
#' FSobj$writeLM(file, direct = "./", LG = NULL, what = NULL, inferGeno = TRUE)
#' ICobj$writeLM(file, direct = "./", LG = NULL, what = NULL, inferGeno = TRUE)
#' 
#' @name $writeLM
#' @param file Character value giving the name of the file to write to.
#' @param direct Character value for the directory to write the file to (relative to the current working working directory)
#' @param LG Integer vector giving the indices of the linkage groups to write. If \code{NULL}, the results for all the 
#' linkage groups are returned.
#' @param what Character vector specifying whether the combined linkage groups \code{"LG-comb"} are to be 
#' returned or the linkage mapping results \code{"map"}.
#' @param inferGeno Logical value indicating whether the paternal and maternal haplotypes of the progeny is to be returned.
#' @return 
#' The function returns a csv file with the rows representing the SNPs and the columns containing the linkage group and linkage map 
#' information.
#' The columns of the file are:
#' \describe{
#' \item{LG}{The linkage number as given in the output from the FS object.}
#' \item{LG_POS}{The position of the SNPs across the linkage group. Starts at 1 and
#' increases by one for each SNP in the same linkage group.}
#' \item{CHROM}{The chromosome name the SNP was located on in the original genomic assembly.}
#' \item{POS}{The position (in base pairs) of the SNP on the original genomic assembly.}
#' \item{TYPE}{The segregation type of the SNP. Either BI (both-informative), 
#' MI (maternal-informative) or PI (paternal-informative).}
#' \item{RF_PAT}{The sum of the recombination fraction estimates for the paternal meioses. Starts at zero
#' for the first SNP in the linkage groups and increases across the linkage group.}
#' \item{RF_MAT}{The sum of the recombination fraction estimates for the maternal meioses. Starts at zero
#' for the first SNP in the linkage groups and increases across the linkage group.}
#' \item{ERR}{The sequencing error estimate for the SNP.}
#' \item{MEAN_DEPTH}{The mean read depth of the SNP.}
#' \item{CALLRATE}{The proportion of individuals in which there is at least one read across the SNP.}
#' }
#' If no linkage maps have been computed (using the \code{\link{$computeMap}}), then the
#' columns RF_PAT, RF_MAT and ERR are not returned. 
#' 
#' If \code{inferGeno=TRUE}, then extra columns are appended to the file containing the paternal and maternal haplotype of the progeny. Each additional column
#' is titled "MAT_[sampleID]" or "PAT_[sampleID]", where "MAT_[sampleID]" denotes the maternal haplotype in individual
#' [sampleID] and "PAT_[sampleID]" is the paternal haplotype in individual [sampleID], where [sampleID] is replaced with the 
#' sampleID of the particular individual. The entries of these extra columns are either \code{A} for the reference allele or \code{B} for the 
#' alternate allele. 
#' 
#' @author Timothy P. Bilton
#' @seealso \code{\link{BC}}, \code{\link{FS}}, \code{\link{IC}}
#' @examples
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(list(sample(c(1,2,4), size=30, replace=TRUE)))
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