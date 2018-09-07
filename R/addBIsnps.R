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
#' Add BI SNPs to existing linkage groups
#' 
#' Method for adding both-informative (BI) SNPs to paternal and maternal linkage groups
#' and to combine maternal and paternal linkage groups.
#' 
#' details here
#' 
#' @usage
#' FSobj$addBIsnps(LODthres = 10, nComp = 30)
#' 
#' @param LODthres An positive numeric value for the LOD threshold used in the 
#' @param nComp An positive integer value 
#' 
#' @author Timothy P. Bilton
#' @name $addBIsnps
NULL