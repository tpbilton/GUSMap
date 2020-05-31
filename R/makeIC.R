##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2020 Timothy P. Bilton <timothy.bilton@agresearch.co.nz>
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
#' Make a intercross (F2) population
#'
#' Create an IC object from an RA object, perform standard filtering and compute statistics specific to intercross populations.
#' 
#' This function converts an RA object into an IC (intercross) object. An IC object is a R6 type object 
#' that contains RA data, various other statistics computed and functions (methods) for analyzing
#' the intercross population and for performing linkage mapping. The statistics computed and data filtering are 
#' specific to intercross populations and sequencing data.
#' 
#' The filtering criteria currently implemented are:
#' \itemize{
#' \item{Minor allele frequency (\code{MAF}): }{SNPs are discarded if their MAF is less than the threshold (default is 0.05)}
#' \item{Proportion of missing data (\code{MISS}): }{SNPs are discarded if the proportion of individuals with no reads (e.g. missing genotype)
#'  is greater than the threshold value (default is 0.5).}
#' \item{Bin size for SNP selection (\code{BIN}): }{SNPs are binned together if the distance (in base pairs) between them is less than the threshold value (default is 100).
#' One SNP is then randomly selected from each bin and retained for final analysis. This filtering is to ensure that there is only one SNP on each sequence read.}
#' \item{Parental read depth (\code{DEPTH}): }{SNPs are discarded if the read depth of either parent is less than the threshold value (default is 5). 
#' This filter is to remove SNPs where the parental information is insufficient to infer segregation type accurately.}
#' \item{Segregation test P-value (\code{PVALUE}): }{SNPs are discarded if the p-value from a segregation test is smaller than the threshold (default is 0.01).
#'  This filters out SNPs where the segregation type has been inferred wrong.}
#' }
#' The segregation type of each SNP is inferred based on the genotypes of the parents. The parental genotypes are called homozygous for the 
#' reference allele if there is only reference reads seen, heterozygous if at least one read for the reference and alternate allele are seen,
#' and homozygous for the alternate allele if only reads for the alternate allele are seen. as a result, the parental genotype may be incorrectly inferred
#' if the read depth is too low (e.g., homozeygous genotype is called heterozygous) and hence why the \code{DEPTH} filter is implemented.
#' The segregation test performed for the \code{PVALUE} filter is described in the supplementary methods 
#' of the publication by \insertCite{bilton2018genetics1;textual}{GUSMap} (Section 4 of File S1).
#' 
#' The pedigree file must be a csv file containing the five columns:
#' \itemize{
#' \item SampleID: A unique character string of the sample ID. These correspond to those found in the VCF file
#' \item IndividualID: A character giving the ID number of the individual for which the sample corresponds to.
#' Note that some samples can be from the same individual. 
#' \item Mother: The ID of the mother as given in the IndividualID. Note, if the mother is unknown then this should be left blank.
#' \item Father: The ID of the father as given in the IndividualID. Note, if the father is unknown then this should be left blank.
#' \item Family: The name of the Family for a group of progeny with the same parents. Note that this is not necessary (it works all the full-sib families) but if
#' given must be the same for all the progeny.
#' }
#' Grandparents can also be supplied but are only used to infer parental genotypes 
#' when the associated read depth is greater than or equal the threshold \code{DEPTH}.
#' 
#' The \code{family} argument allows the user to specify the family to be used in the creation of
#' the full-sib population. Note that this argument using the "Family" column in the pedigree file and
#' so pedigree file needs to be set-up correctly.
#' 
#' Note: Only a single intercross family can be processed at present. There are future plans to extend this out to include
#' multiple families.
#' 
#' @param RAobj Object of class RA created via the \code{\link[GUSbase]{readRA}} function.
#' @param pedfile Character string giving the file name (relative to the current directory) of the pedigree file.
#' @param family Vector of character strings giving the families to retain in the IC object. This allows a pedigree file with more than one family to be supplied.
#' @param filter Named list of thresholds for various filtering criteria.
#' See below for details.
#'  
#' @return 
#' An R6 object of class \code{\link{IC}}
#' @author Timothy P. Bilton
#' @examples 
#' ## extract filename for Manuka dataset in GUSMap package
#' vcffile <- Manuka11()
#' 
#' ## Convert VCF to RA format
#' rafile <- VCFtoRA(vcffile$vcf)
#' 
#' ## read in the RA data
#' mkdata <- readRA(rafile)
#' 
#' ## Create the IC population
#' makeIC(mkdata, pedfile=vcffile$ped)
#' @references
#' \insertRef{bilton2018genetics1}{GUSMap}
#' @export

### Make a full-sib family population
makeIC <- function(RAobj, samID=NULL,
                   filter=list(MAF=0.05, MISS=0.2, BIN=100, DEPTH=5, PVALUE=0.01, MAXDEPTH=500)){
  #inferSNPs = FALSE, perInfFam=1){
  inferSNPs = FALSE; perInfFam=1 # some variables for multiple familes needed for later
  ## Do some checks
  if(!all(class(RAobj) %in% c("RA","R6")))
    stop("First argument supplied is not of class 'R6' and 'RA'")
  if(is.null(filter$MAF) || GUSbase::checkVector(filter$MAF, type="pos_numeric", minv=0, maxv=1) || length(filter$MAF) != 1){
    message("Minor allele filter has not be specified or is invalid. Setting to 5%:")
    filter$MAF <- 0.05
  }
  if(perInfFam <= 0.5 || perInfFam > 1)
    stop("The of the percentage of families which are informative for each SNP must greater than 50% and leass than equal to 100%.")
  if(is.null(filter$MISS) || GUSbase::checkVector(filter$MISS, type="pos_numeric", minv=0, maxv=1) || length(filter$MISS) != 1){
    message("Proportion of missing data filter has not be specified or is invalid. Setting to 20%:")
    filter$MISS <- 0.2
  }
  if(is.null(filter$BIN) || GUSbase::checkVector(filter$BIN, type="pos_numeric", minv=0) || length(filter$BIN) != 1){
    message("Minimum distance between adjacent SNPs is not specified or is invalid. Setting to 0:")
    filter$BIN <- 0 
  }
  if(is.null(filter$DEPTH) || GUSbase::checkVector(filter$DEPTH, type="pos_numeric", minv=0) || length(filter$DEPTH) != 1){
    message("Minimum depth on the parental genotypes filter has not be specified or is invalid. Setting to a depth of 5")
    filter$DEPTH <- 5
  }
  if(is.null(filter$PVALUE) || GUSbase::checkVector(filter$PVALUE, type="pos_numeric", minv=0, maxv = 1) || length(filter$PVALUE) != 1){
    message("P-value for segregation test is not specified or invalid. Setting a P-value of 0.01:")
    filter$PVALUE <- 0.05
  }
  if(is.null(filter$MAXDEPTH) || GUSbase::checkVector(filter$MAXDEPTH, type="pos_numeric", minv=0) || length(filter$MAXDEPTH) != 1){
    message("Maximum average SNP read depth is not specified or invalid. Setting to 500:")
    filter$MAXDEPTH <- 500
  }
  
  ## initalize the UR object
  ICobj <- IC$new(RAobj)

  ## Extract the private variables we want
  indID <- ICobj$.__enclos_env__$private$indID
  nSnps <- ICobj$.__enclos_env__$private$nSnps
  
  if(is.null(samID)) samID = indID
  else {
    samID = unique(samID)
    if(!is.vector(samID) || !is.character(samID) || any(is.na(samID)) || any(!(samID %in% indID)))
      stop( paste0("Sample IDs ", paste0(which(!(samID %in% indID)), collapse = ","), " are not present in the dataset"))
  }
  indx_ind = which(indID %in% samID)
  genon = ICobj$.__enclos_env__$private$genon[indx_ind,]
  ref   = ICobj$.__enclos_env__$private$ref[indx_ind,]
  alt   = ICobj$.__enclos_env__$private$alt[indx_ind,]
  nInd  = nrow(ref)
  indID = indID[indx_ind]
  
  ## Define variables that will be used.
  cat("-------------\n")
  cat("Processing Data.\n\n")
  
  cat("Filtering criteria for removing SNPs:\n")
  cat("Minor allele frequency (MAF) < ", filter$MAF,"\n", sep="")
  cat("Percentage of missing genotypes > ", filter$MISS*100,"%\n", sep="")
  cat("Distance for binning SNPs <= ", filter$BIN," bp\n", sep="")
  cat("Read depth associated with at least one parental genotype <= ", filter$DEPTH,"\n", sep="")
  cat("P-value for segregation test < ", filter$PVALUE,"%\n", sep="")
  cat("Average SNP depth is > ", filter$MAXDEPTH,"\n\n", sep="")

  
  ## check the pedigree
  # if(!is.null(pedfile)){
  #   ped = utils::read.csv(pedfile, stringsAsFactors=F)
  #   if(GUSbase::checkVector(family, type="one_logical") || !(family %in% ped$Family))
  #     stop("Family name is not supplied")
  #   indx_ind = which(ped[,5] == family)
  #   genon = ICobj$.__enclos_env__$private$genon[indx_ind,]
  #   ref   = ICobj$.__enclos_env__$private$ref[indx_ind,]
  #   alt   = ICobj$.__enclos_env__$private$alt[indx_ind,]
  # } else{
  #   genon = ICobj$.__enclos_env__$private$genon
  #   ref   = ICobj$.__enclos_env__$private$ref
  #   alt   = ICobj$.__enclos_env__$private$alt
  # }
    
  noFam <- as.integer(1); fam=1
  config_all <- nInd_all <- indx <- indID_all <- vector(mode = "list", length = noFam)
  genon_all <- ref_all <- alt_all <- vector(mode="list", length=noFam)
  
  ## Run the filtering of the progeny SNPs
  MAF <- colMeans(genon, na.rm=T)/2
  MAF <- pmin(MAF,1-MAF)
  miss <- apply(genon,2, function(x) sum(is.na(x))/length(x))
  maxdepth <- colMeans(ref + alt)
  indx_temp <- rep(TRUE, nSnps)
  indx_temp[which(MAF < filter$MAF | miss > filter$MISS | maxdepth > filter$MAXDEPTH)] <- FALSE
  
  ## fix max depth to be at 500 total reads
  ## Run into numeric issues otherwise
  high_depth <- which((ref+alt) > 500)
  ref_temp <- ref[high_depth]
  alt_temp <- alt[high_depth]
  ref[high_depth] <- as.integer(ref_temp/(ref_temp+alt_temp)*500)
  alt[high_depth] <- as.integer(alt_temp/(ref_temp+alt_temp)*500)
  
  ## assume segregation type is BI
  config <- rep(1,nSnps)
  
  #### Segregation test to determine if the SNPs have been miss-classified
  seg_Dis <- sapply(1:nSnps,function(x){
    if(is.na(config[x]))
      return(NA)
    else{
      d = ref[,x] + alt[,x]
      g = genon[,x]
      K = sum(1/2^(d[which(d != 0)])*0.5)/sum(d != 0)
      nAA = sum(g==2, na.rm=T)
      nAB = sum(g==1, na.rm=T)
      nBB = sum(g==0, na.rm=T)
      ## check that there are sufficient data to perform the chisq test
      if(sum(nAA+nAB+nBB)/length(g) <= (1-filter$MISS))
        return(NA)
      else if(config[x] == 1){
        exp_prob <- c(0.25 + K,0.5 - 2*K, 0.25 + K)
        ctest <- suppressWarnings(stats::chisq.test(c(nBB,nAB,nAA), p = exp_prob))
        return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
      }
      else return(NA)
    }
  },simplify = T)
  config[which(seg_Dis)] <- NA
  
  chrom <- ICobj$.__enclos_env__$private$chrom
  chrom[!indx_temp | is.na(config)] <- NA
  pos <- ICobj$.__enclos_env__$private$pos
  pos[!indx_temp | is.na(config)] <- NA
  ## Extract one SNP from each read.
  set.seed(36475)
  if(filter$BIN > 0){
    oneSNP <- rep(FALSE,nSnps)
    oneSNP[unlist(sapply(unique(chrom[which(!is.na(chrom))]), function(x){
      ind <- which(chrom == x)
      g1_diff <- diff(pos[ind])
      SNP_bin <- c(0,cumsum(g1_diff > filter$BIN)) + 1
      set.seed(58473+as.numeric(which(x==chrom))[1])
      keepPos <- sapply(unique(SNP_bin), function(y) {
        ind2 <- which(SNP_bin == y)
        if(length(ind2) > 1)
          return(sample(ind2,size=1))
        else if(length(ind2) == 1)
          return(ind2)
      })
      return(ind[keepPos])
    },USE.NAMES = F ))] <- TRUE
  }
  else 
    oneSNP <- rep(TRUE,nSnps)
  
  indx_temp[which(indx_temp & (!oneSNP | is.na(config)))] <- FALSE
  indx[[fam]] <- indx_temp
  
  ## Determine the segregation groups
  config_all[[fam]] <- config
  nInd_all[[fam]] <- nInd
  genon_all[[fam]] <- genon
  ref_all[[fam]] <- ref
  alt_all[[fam]] <- alt
  
  if(noFam == 1){
    fam = 1
    
    ## Find all the SNPs to keep and subset the global variables
    indx_all <- indx[[fam]]
    indID_all[[fam]] <- indID

    genon_all[[fam]]     <- genon[,indx_all]
    ref_all[[fam]]       <- ref[,indx_all]
    alt_all[[fam]]       <- alt[,indx_all]
    chrom_all            <- ICobj$.__enclos_env__$private$chrom[indx_all]
    pos_all              <- ICobj$.__enclos_env__$private$pos[indx_all]
    SNP_Names            <- ICobj$.__enclos_env__$private$SNP_Names[indx_all]
    
    #group <- group_infer <- vector(mode="list", length=noFam)
    group <- list()
    group$BI <- which(config_all[[fam]][indx_all] == 1)

    config_all[[fam]] <- config_all[[fam]][indx_all]

    cat("-------------\n")
    cat("Summary:\n\n")
    cat("Number of SNPs remaining after filtering:",sum(indx_all),"\n")
    cat("Number of progeny:", nInd_all[[fam]],"\n")
    
    ## Create the summary info:
    temp <- ref_all[[1]] + alt_all[[1]]
    summaryInfo <- list()
    summaryInfo$data <- paste0(c(
      "Single Family Linkage analysis:\n\n",
      "Data Summary:\n",
      "Data file:\t", RAobj$.__enclos_env__$private$infilename,"\n",
      "Mean Depth:\t", round(mean(temp),4),"\n",
      "Mean Call Rate:\t",round(sum(temp!=0)/length(temp),4),"\n",
      "Number of ...\n",
      "  Progeny:\t",unlist(nInd),"\n",
      "  SNPs:\t\t",length(unlist(group)),"\n\n"))
  }
  else{
    stop("yet to be implemented")
  }
  
  ## Update the R6 object and return it
  ICobj$.__enclos_env__$private$updatePrivate(list(
    genon = genon_all, ref = ref_all, alt = alt_all, chrom = chrom_all, pos = pos_all,
    group = group, config = config_all, nInd = nInd_all, nSnps = sum(indx_all), 
    noFam = noFam, indID = indID_all, SNP_Names = SNP_Names,
    masked=rep(FALSE,sum(indx_all)), famInfo=NULL, summaryInfo=summaryInfo)
  )
  
  return(ICobj)
}
