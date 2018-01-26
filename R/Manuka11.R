##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping
# Copyright 2017 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
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
### R Function for reading in the Manuka data of chromosome 11.
### Author: Timothy Bilton
### Date: 26/01/18

## Wrapper function for reading in the Manuka data set for chromosome 11


#' Manuka chromosome 11 SNPs
#' 
#' Function for extracting the path to the (VCF) file of the Manuka data used in the publication by Bilton et
#' al. (2017).
#' 
#' The data consists of 680 SNPs, genotyped using genotyping-by-sequencing
#' methods. The data is in VCF format (see this \href{https://samtools.github.io/hts-specs/VCFv4.2.pdf}{page} 
#' for specification of VCF format).
#' 
#' @return Function outputs a character string of the complete path to the data set in the package.
#' @author Timothy P. Bilton
#' @references Bilton, T.P., Schofield, M.R., Black, M.A., Chagné, D., Wilcox,
#' P.L., Dodds K.G. (2017). Accounting for errors in low coverage high-throughput
#' sequencing data when constructing genetic maps using biparental outcrossed
#' populations. Unpublished Manuscript.
#' @examples
#' 
#' ## extract the name of the vcf file 
#' Manuka11()
 

#' 
#' @export Manuka11
Manuka11 <- function(){
  return(system.file("extdata", "Manuka11_chr11.vcf", package="GUSMap"))
}

