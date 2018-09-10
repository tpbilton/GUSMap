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
#' FS method: Print summary information
#' 
#' Method for printing summary information of an FS object
#' 
#' There are four types of information availabe to be printed for an object (as specified by the \code{what} argument).
#' \describe{
#' \item{\code{"data"}}{Summary information of the data as process by the \code{\link{makeFS}} function.}
#' \item{\code{"LG-pts"}}{Summary information of the pesudo-testcross linkage groups produced by the \code{\link{$createLG}} function.}
#' \item{\code{"LG-comb"}}{Summary information of the combined linkage groups produced by the \code{\link{$addBIsnps}} function.}
#' \item{\code{"map"}}{Summary information of the linkage maps produced by the \code{\link{$computeMap}} function.}
#' }
#' 
#' @usage
#' FSobj$print(what = NULL, ...)
#' 
#' @param what Character vector giving the type of information to be printed. IF \code{NULL}, then all the 
#' information is printed.
#' 
#' @name $print
#' @author Timothy P. Bilton
#' @seealso \code{\link{FS}}
#' @examples 
#' ## simulate some sequencing data
#' set.seed(6745)
#' config <- list(replicate(2,sample(c(1,2,4), size=30, replace=T), simplify=FALSE))
#' F1data <- simFS(0.01, config=config, meanDepth=10, nInd=50)
#' # print summary information of data
#' F1data$print(what = "data")
#' 
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt()
#' ## create paternal and maternal linkage groups
#' F1data$createLG()
#' ## Add BI SNPs
#' F1data$addBIsnps()
#' 
#' # print information of the combined linkage groups
#' F1data$print(what = "LG-comb")
#' #print both the data infromation and separate linkage groups
#' F1data$print(what = c("data", "LG-pts"))
NULL