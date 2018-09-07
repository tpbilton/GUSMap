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
#' Order linkage groups
#' 
#' Method for ordering linakge groups in GUSMap
#' 
#' Linkage groups created using the \code{\link{$addBIsnps}} are ordered using the multidimensional scaling (MDS)
#' approach described by \insertCite{preedy2016tag;textual}{GUSMap}. 
#' 
#' 
#' @usage
#' FSobj$orderLG(chrom = NULL, mapfun = "haldane", weight="LOD2", ndim=30)
#' 
#' @param chrom An integer vector of the indices for the linkage group(s) that are to be ordered.
#' @param mapfun A character value for the mapping function to be used 
#' @param weight A character value for the weight function 
#' @param ndim A integer value for the number of dimensions to use in the MDS 
#' 
#' @name $orderLG
#' @references 
#' \insertRef{preedy2016tag}{GUSMap} 
NULL