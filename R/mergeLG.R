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
#' Merge Linkage Groups
#' 
#' Method for merging multiple linkage group from the list of linkage groups formed
#' in an FS object.
#' 
#' In a linkage analysis in GUSMap, there may be need to remove linkage groups from the 
#' list of created linkage groups. The indices of the linkage groups corresponds to the number
#' given in the output for the FS object. There needs to be at leads two linakge groups 
#' specified for merging and more than two linkage groups can be merged at once.
#' 
#' @usage
#' FSobj$mergeLG(LG)
#' 
#' @param LG An integer vector specifying the indices of the linkage groups to be merged. 
#' @name $mergeLG
#' @author Timothy P. Bilton
NULL