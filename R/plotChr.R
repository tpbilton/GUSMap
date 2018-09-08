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
#' FS method: Heatmap of 2-point rf and LOD estimates from chromosome order
#' 
#' Method for plotting 2-point recombination fraction and LOD score estimates when the SNPs are ordered according to the genomic assembly.
#'
#' This function plots the heatmap of the matrix of 2-point recombination fraction estimates (or 2-point LOD scores)
#' as computed from the \code{\link{$rf_2pt}} function when the SNPs are ordered to the genomic assembly (as given in the VCF file). 
#' 
#' 
#' Note: 
#' 
#' @usage
#' FSobj$plotChr(parent = "maternal", mat="rf", filename=NULL, chrS=2, lmai=2)
#' 
#' @param parent Character value specifying whether the SNPs segreagting in the maternal parent should be 
#' plotted (\code{"maternal"}), or whether whether the SNPs segreagting in the paternal parent should be 
#' plotted (\code{"paternal"}), or whether all the SNPs should be plotted (\code{"both"}).
#' @param mat Charater value for the matrix to be plotted. \code{"rf"} plots the matrix of 2-point 
#' recombination fractions while \code{"LOD"} plots the matrix of 2-point LOD scores.
#' @param filename Character giving the name of the file to save the plot to. If \code{NULL}, the plot is given 
#' in the graphics window and not saved to a file.
#' @param chrS Numeric value. Controls the size of the chromosome names on the left side of the plot.
#' @param lmai Numeric value. Controls the amount of white space for the chromosome names on
#' the left side of the plot.
#' 
#' @name $plotChr
#' @author Timothy P. Bilton
#' @seealso \code{\link{FS}}
#' @examples 
#' ## Simulate some sequencing data
#' set.seed(6745)
#' config <- list(list(sample(c(1,2,4), size=10, replace=T)), list(sample(c(1,2,4), size=10, replace=T)))
#' F1data <- simFS(0.01, config=config, meanDepth=10, 
#' ## Compute 2-point recombination fractions
#' F1data$rf_2pt()
#' 
#' ## Plot the heatmap
#' F1data$plotChr()
#' 
#' @aliases NULL
