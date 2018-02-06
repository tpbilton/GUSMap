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
#### Functions for converting VCF files to RA data
#### Author: Timothy P. Bilton
#### Adapted from a python script written by Rudiger Brauning and Rachael Ashby 
#### Date: 06/02/18

#' Convert VCF file into RA (Reference/Alternative) file.
#'
#' Function for converting a VCF file into RA format.
#' 
#' The VCF files must contain some information regarding allelic depth. Currently, the function can use one of the following
#' fields (or group of fields) in a VCF file:
#' \itemize{
#' \item AD field
#' \item AO and RO fields
#' \item DP4 field
#' }
#' Information regarding VCF files and their format can be found at the samtools \href{https://samtools.github.io/hts-specs/VCFv4.3.pdf}{github} page.
#' 
#' RA format is a tab-delimited with columns, CHROM, POS, SAMPLES
#' where SAMPLES consists of sampleIDs, which typically consist of a colon-delimited sampleID, flowcellID, lane, seqlibID.
#' e.g.,
#' \tabular{llll}{
#' CHROM \tab  POS  \tab   999220:C4TWKACXX:7:56 \tab  999204:C4TWKACXX:7:56 \cr
#' 1     \tab  415  \tab   5,0                   \tab  0,3                   \cr
#' 1     \tab  443  \tab   1,0                   \tab  4,4                   \cr
#' 1     \tab  448  \tab   0,0                   \tab  0,2
#' }
#' Note: Indels are removed, multiple alternative alleles are removed and ./. is translated into 0,0.
#' 
#' The format of the pedigree files is a csv file with the following columns.
#' \itemize{
#' \item SampleID: A unique character string of the sample ID. These must correspond to those found in the VCF file.
#' \item IndividualID: A character giving the ID number of the individual for which the sample corresponds to.
#' Note that some samples can be from the same individual. 
#' \item Mother: The ID of the mother as given in the IndividualID. Note, if the mother is unknown then this should be left blank.
#' \item Father: The ID of the father as given in the IndividualID. Note, if the father is unknown then this should be left blank.
#' \item Family: The name of the Family for a group of progeny with the same parents. Note that this is not necessary but if
#' given must be the same for all the progeny.
#' }
#' 
#' @param infilename String giving the filename of the VCF file to be converted to RA format
#' @param direct String of the directory (or relative to the working direct) where the RA file is to be written.
#' @param makePed A logical value. If TRUE, a pedigree file is initialized.
#' @return A string of the complete file path and name of the RA file created from the function.
#' In addition to creating a RA file, a pedigree file is also initialized in the same folder as the RA file if
#' specified and the named pedigree does not already exist.
#' @author Timothy P. Bilton. Adapted from a Python script written by Rudiger Brauning and Rachael Ashby.
#' @seealso \code{\link{readRA}}
#' @examples
#' MKfile <- Manuka11()
#' RAfile <- VCFtoRA(MKfile$vcf, makePed=F)
#' @export VCFtoRA

VCFtoRA <- function(infilename, direct="./", makePed=T){
  
  ## Do some checks
  if(!is.character(infilename) || length(infilename) !=1)
    stop("The input file name is not a string of length 1.")
  if(!file.exists(infilename))
    stop("Input file does not exist. Check your wording or the file path.")
  if(!is.character(direct) || length(direct) != 1)
    stop("Invalid input for the path to the directory where the RA file is to be written.")

  cat("Processing VCF file: Converting to RA format.\n\n")
    
  outfilename <- paste0(tail(strsplit(infilename,split=.Platform$file.sep)[[1]],1),".ra.tab")
  outpath <- dts(normalizePath(direct, winslash=.Platform$file.sep, mustWork=T))
  
  headerlist = c('CHROM', 'POS')
  
  ## Read in the lines of the file
  Lines <- readLines(infilename)
  
  ## entries for empty genotypes
  empty_genotypes <- c("./.",".,.",".",".|.")
  
  start = which(unlist(lapply(Lines, function(x) substr(x,1,2) == "##")))  ## starting position
  start = max(start,0) + 1
  
  ## create the file
  outfile = file.path(outpath,outfilename)
  file.create(outfile,showWarnings = F)

  ## write the heading line
  line = strsplit(Lines[start], split="\t")[[1]]
  headerlist = c(headerlist, line[10:length(line)])
  ## Write the header to the RA file
  newLines <- list(paste(headerlist, collapse="\t"))
  cat("Found",length(headerlist),"samples\n")
  
  ## Now write the SNPs
  for(i in (start+1):length(Lines)){
    line = trimws(Lines[[i]])

    line = strsplit(line, split="\t")[[1]]
    chrom = line[1]
    pos = line[2]
    ref = line[4]
    alt = line[5]
    ## Filter out indels and multiple alternative alleles
    if (ref == "." || alt == "." || length(alt) > 1)
      next
    else{
      ## Check that there is allelic depth in the VCF file
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
      for(j in line[10:length(line)]){
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
      newLines[[i-(start)+1]] <- paste0(c(chrom,pos,newline), collapse = "\t")
    }
  }
  ## open the connection to the file
  con <- file(outfile)
  writeLines(unlist(newLines), con=con)
  close(con)
  ## output the information
  cat(length(Lines)-start,"SNPs written\n\n")
  cat("Name of RA file:    ",outfilename,"\n")
  cat("Location of RA file: ",outpath,"/\n\n",sep="")
  ## Initialize the pedigree file
  if(makePed){
    pedfile <- paste0(strsplit(paste0(tail(strsplit(infilename,split=.Platform$file.sep)[[1]],1)),split="\\.")[[1]][1],"_ped.csv")
    pedpath <- file.path(outpath,pedfile)
    if(!file.exist(pedpath)){
      cat("A pedigree file has been initialized.\n")
      cat("Name of pedigree file:     ",pedfile,"\n")
      cat("Location of pedigree file: ",outpath,"/\n\n",sep="")
      nSamp <- length(headerlist) - 2
      write.csv(cbind("SampleID"=c(headerlist[-c(1,2)]),"IndividualID"=rep("",nSamp),"Mother"=rep("",nSamp),
                  "Father"=rep("",nSamp), "Family"=rep("",nSamp)), file = pedpath, row.names=F)
    }
  }  
  return(invisible(outfile))
}

#### Some functions from the kutils package for removing trailing spaces for filenames.
dts <- function (name) 
  gsub("/$", "", dms(name))

dms <- function(name)
  gsub("(/)\\1+", "/", name)

