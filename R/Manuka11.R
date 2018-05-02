##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping
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
### R Function for reading in the Manuka data of chromosome 11.
### Author: Timothy Bilton
### Date: 6/02/18

## Wrapper function for reading in the Manuka data set for chromosome 11
#' Manuka chromosome 11 SNPs and pedigree file
#' 
#' Function for extracting the path to the (VCF) file of the Manuka data used in the publication by Bilton et
#' al. (2018) and a pedigree file that goes with the data.
#' 
#' The data consists of 680 SNPs, genotyped using genotyping-by-sequencing
#' methods. The data is in VCF format (see this \href{https://samtools.github.io/hts-specs/VCFv4.2.pdf}{page} 
#' for specification of VCF format). The pedigree file contains five columns, namely
#' \itemize{
#' \item SampleID: A unique character string of the sample ID. These correspond to those found in the VCF file
#' \item IndividualID: A character giving the ID number of the individual for which the sample corresponds to.
#' Note that some samples can be from the same individual. 
#' \item Mother: The ID of the mother as given in the IndividualID. Note, if the mother is unknown then this should be left blank.
#' \item Father: The ID of the father as given in the IndividualID. Note, if the father is unknown then this should be left blank.
#' \item Family: The name of the Family for a group of progeny with the same parents. Note that this is not necessary but if
#' given must be the same for all the progeny.
#' }
#' 
#' @return Function outputs a character string of the complete path to the manuka data set contained in the package
#' and a pedigree file that goes with the data.
#' @author Timothy P. Bilton
#' @references Bilton, T.P., Schofield, M.R., Black, M.A., Chagné, D., Wilcox,
#' P.L., & Dodds, K.G. (2018). Accounting for errors in low coverage high-throughput
#' sequencing data when constructing genetic maps using biparental outcrossed
#' populations. \emph{Genetics}, \emph{209}(1), 65--76. \url{http://www.genetics.org/content/209/1/65}
#' @examples
#' ## extract the name of the vcf fileand the pedigree file
#' Manuka11()
#' 
#' @export Manuka11

Manuka11 <- function(){
  return(list(vcf=system.file("extdata", "Manuka_chr11.vcf", package="GUSMap"),
              ped=system.file("extdata", "Manuka_chr11_ped.csv", package="GUSMap")))
}

