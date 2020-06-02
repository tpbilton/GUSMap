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
#' BC, FS and IC method: Plot linkage maps
#' 
#' Method for plotting linkage maps in GUSMap
#' 
#' Function displays linkage maps computed using the \code{\link{$computeMap}} function.
#' 
#' For a list of mapping functions available, see \code{\link{mfun}}. 
#'  
#' @usage
#' BCobj$plotLM(LG = NULL, fun="haldane", col="black")
#' FSobj$plotLM(LG = NULL, fun="haldane", col="black")
#' ICobj$plotLM(LG = NULL, fun="haldane", col="black")
#' 
#' @param LG An integer vector specifying the indices of the linkage groups for which their linkage maps are to be plotted. 
#' @param fun Character value giving the mapping function to use which is passed to \code{\link{mfun}}.
#' @param col Vector of character values giving the colors of the linkage maps.
#' 
#' @name $plotLM
#' @author Timothy P. Bilton
#' @seealso \code{\link{BC}}, \code{\link{FS}}, \code{\link{IC}}, \code{\link{mfun}}
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
#' ## Compute Map
#' F1data$computeMap()
#' 
#' ## Produce linkage maps
#' F1data$plotLM()
NULL

