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
#### Functions for converting phase information to OPGPs
#### Author: Timothy P. Bilton
#### Adapted from a python script written by Rudiger Brauning and Rachael Ashby 
#### Date: 15/01/18

#' Convert VCF file into RA (ReferenceAlternative) file.
#'
#' Function for converting a VCF file into RA format.
#'
#' RA format is a tab-delimited with columns, CHROM, POS, SAMPLES
#' where SAMPLES consists of colon-delimited sampleID, flowcellID, lane, seqlibID.
#' e.g.,
#' CHROM   POS     999220:C4TWKACXX:7:56   999204:C4TWKACXX:7:56
#' 1       415     0,0                     0,0
#' 1       443     1,0                     9,0
#' 1       448     0,0                     0,0
#' 
#' Note: Indels are removed, multiple alternative alleles are removed and ./. is translated into 0,0.
#' 
#' @param infilename String giving the filename of the VCF file to be converted to RA format
#' @return Nothing is returned. Instead a new file is written to the same directory as where the 
#' VCF file is found. Note: the name of the RA file outputted is '\code{infilename}.ra.tab'
#' @author Timothy P. Bilton. Adapted from a Python script written by Rudiger Brauning and Rachael Ashby.
#' @export VCFtoRA

VCFtoRA <- function(infilename){
  
  if(is.string(infilename) || length(infilename))
    stop("The input file name is not a string of length 1.")

  outfilename <- paste0(infilename,".ra.tab")
  
  headerlist = c('CHROM', 'POS')
  
  ## Read in the lines of the file
  Lines <- readLines(infilename)

  ## entries for empty genotypes
  empty_genotypes <- c("./.",".,.",".",".|.")
  
  ## start of the SNPs 
  snp_start = 0
  
  cat("Processing VCF file: Converting to RA format.\n\n")
  
  for(i in 1:length(Lines)){
    line = Lines[[i]]
    line = trimws(line)
    
    ## If comments, ignore
    if(substr(line,1,2) == "##")
      next
    ## if the headings nad sample IDS
    else if(substr(line,1,6) == "#CHROM"){
      ## extract the sample IDS
      line = strsplit(line, split="\t")[[1]]
      headerlist = c(headerlist, line[10:length(line)])
      ## Write the header to the RA file
      write(paste(paste(headerlist, collapse="\t"),sep="\n"), file=outfilename, ncolumns = 1)
      cat("found",length(headerlist),"samples\n")
      snp_start = i + 1
    }
    ## If the genotype data
    else{
      line = strsplit(line, split="\t")[[1]]
      chrom = line[1]
      pos = line[2]
      ref = line[4]
      alt = line[5]
      ## Filter out indels and multiple alternative alleles
      if (ref == "." || alt == "." || length(alt) > 1)
        next
      else{
        format = strsplit(line[9], split= ":")[[1]]
        if("AD" %in% format)
          ad_pos = which(format == "AD")
        else if(all(c("RO","AO") %in% format)){
          ro_pos = which(format == "RO")
          ao_pos = which(format == "AO")
        }
        else if("DP4" %in% format)
          dp4_pos = which(format == "DP4")
        else
          stop("We can't use this vcf file. AD (alleleic depth or RO (Reference allele observation count) and AO (Alternate allele observation count) information is needed.\n")
        ## Extract the alleles depths
        newline = c()
        #for(j in line[10:length(line)]){
        for(j in line[10:110]){
          if (j %in% empty_genotypes)
            newline = c(newline,"0,0")
          else{
            j = strsplit(j, split = ":")[[1]]
            if( length(ad_pos) > 0 ){              #If ad_pos is not null, i.e., there is a value for it, append it to outlist
              if( j[ad_pos] %in% empty_genotypes ) #gatk vcf4.2 will fill out genotype fields, even for uncovered data
                newline = c(newline,"0,0")
              else
                newline = c(newline, j[ad_pos])
            }
            else if(length(ro_pos) > 0 && length(ao_pos) > 0){ #ELSE IF ro_pos and ao_pos are not null or are equal to 0
              ad = paste0(j[ro_pos], ",", j[ao_pos])
              newline = c(newline, ad)
            }
            else if( dp4_pos ){
              counts = strsplit(j[dp4_pos], split=",")
              allele1 = as.integer(counts[1]) + as.integer(counts[2])
              allele2 = as.integer(counts[3]) + as.integer(counts[4])
              ad = paste0(allele1, ",", allele2)
              newline = c(newline, ad)
            }
            else ##Should never really get here, but if AD, AO and RO are all null, it will break the script
              stop("Can't Find either Allele Depth (AD) or RO (Reference allele observation count) and AO (Alternate allele observation count) at this position.\n")
          }
        }
        write(paste(newline, collapse="\t"), file = outfilename, append = T)
      }
    }
  }
  cat(length(Lines)-snp_start,"SNPs written\n\n")
}

