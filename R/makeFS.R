##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2023 Timothy P. Bilton 
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
#' Make a full-sib family (FS) population
#'
#' Create an FS object from an RA object, perform standard filtering and compute statistics specific to full-sib family populations.
#' 
#' This function converts an RA object into an FS (full-sib) object. An FS object is a R6 type object 
#' that contains RA data, various other statistics computed and functions (methods) for analyzing and performing linkage mapping
#' on full-sib families. The statistics computed and data filtering are specific to full-sib family populations and sequencing data.
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
#' \item{Maximum averge SNP depth (\code{MAXDEPTH}}{SNPs with an average read depth above the threshold value are discarded.}
#' }
#' The segregation type of each SNP is inferred based on the genotypes of the parents. The parental genotypes are called homozygous for the 
#' reference allele if there is only reference reads seen, heterozygous if at least one read for the reference and alternate allele are seen,
#' and homozygous for the alternate allele if only reads for the alternate allele are seen. As a result, the parental genotype may be incorrectly inferred
#' if the read depth is too low (e.g., homozeygous genotype is called heterozygous) and hence why the \code{DEPTH} filter is implemented.
#' The segregation test performed for the \code{PVALUE} filter is described in the supplementary methods 
#' of the publication by \insertCite{bilton2018genetics1;textual}{GUSMap} (Section 4 of File S1).
#' 
#' If the argument \code{inferSNPs} is \code{TRUE}, an attempt to infer the segregation type of SNPs where the segregation type could 
#' not be determined from the parental genotypes is made using the progeny data only. Note that using this approach, MI SNPs can not be 
#' distinguished from PI SNPs (since we only know that one parent is heterozygous and one parent is homozygous but we don't know which is which)
#' and so we collectively refer to the MI and PI SNPs inferred using this approach as semi-informative (SI) SNPs.
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
#' The \code{keepSNPs} will over-ride any filtering criteria specified in \code{filter} but does require that the segregation type can be inferred. 
#' That is, a SNP will be retained in the dataset if its indices is in the \code{keepSNPs} vector even though it may fail one or more of the filtering 
#' criteria as defined in \code{filter}, except in the case where the segregation type cannot be inferred.
#' 
#' Note: Only a single full-sib family can be processed at present. There are future plans to extend this out to include
#' multiple families.
#' 
#' @param RAobj Object of class RA created via the \code{\link[GUSbase]{readRA}} function.
#' @param pedfile Character string giving the file name (relative to the current directory) of the pedigree file.
#' @param family Vector of character strings giving the families to retain in the FS object. This allows a pedigree file with more than one family to be supplied.
#' @param inferSNPs Logical value indicating whether to infer the segregation type of SNPs using the progeny information only
#' in cases where the segregation type could not be inferred from the parental genotypes.
#' @param keepSNPs Vector of indices of the SNPs in the RA to retain in the full-sib population.
#' @param filter Named list of thresholds for various filtering criteria.
#' See below for details.
#'  
#' @return 
#' An R6 object of class \code{\link{FS}}
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
#' ## Create the FS population
#' makeFS(mkdata, pedfile=vcffile$ped)
#' @references
#' \insertRef{bilton2018genetics1}{GUSMap}
#' @export

#### Note: need to add a check that the SNPs are in order (within a chromosome) somewhere:


### Make a full-sib family population
makeFS <- function(RAobj, pedfile, family=NULL, inferSNPs=TRUE, keepSNPs=NULL,
                   filter=list(MAF=0.05, MISS=0.2, BIN=100, DEPTH=5, PVALUE=0.01, MAXDEPTH=500)){
  #inferSNPs = FALSE, perInfFam=1){
  perInfFam=1; MNIF=1 # some variables for multiple families needed for later
  ## Do some checks
  if(!all(class(RAobj) %in% c("RA","R6")))
    stop("First argument supplied is not of class 'R6' and 'RA'")
  if(is.null(filter$MAF) || GUSbase::checkVector(filter$MISS, minv=0, maxv=1, type="pos_numeric", equal=FALSE) || length(filter$MISS) != 1){
    warning("Minor allele filter has not be specified or is invalid. Setting to 5%:")
    filter$MAF <- 0.05
  }
  if(perInfFam <= 0.5 || perInfFam > 1)
    stop("The of the percentage of families which are informative for each SNP must greater than 50% and leass than equal to 100%.")
  if(is.null(filter$MISS) ||  GUSbase::checkVector(filter$MISS, minv=0, maxv=1, type="pos_numeric", equal=FALSE) || length(filter$MISS) != 1){
    warning("Proportion of missing data filter has not be specified or is invalid. Setting to 20%:")
    filter$MISS <- 0.2
  }
  if(is.null(filter$BIN) || GUSbase::checkVector(filter$BIN, minv=0, type="pos_numeric", equal=FALSE) || length(filter$BIN) != 1){
    warning("Minimum distance between adjacent SNPs is not specified or is invalid. Setting to 0:")
    filter$BIN <- 0 
  }
  if(is.null(filter$DEPTH) || GUSbase::checkVector(filter$DEPTH, minv=0, type="pos_numeric", equal=FALSE) || length(filter$DEPTH) != 1){
    warning("Minimum depth on the parental genotypes filter has not be specified or is invalid. Setting to a depth of 5")
    filter$DEPTH <- 5
  }
  if(is.null(filter$PVALUE) || GUSbase::checkVector(filter$PVALUE, minv=0, maxv=1, type="pos_numeric", equal=FALSE) || length(filter$PVALUE) != 1){
    warning("P-value for segregation test is not specified or invalid. Setting a P-value of 0.01:")
    filter$PVALUE <- 0.01
  }
  if(is.null(filter$MAXDEPTH) || GUSbase::checkVector(filter$MAXDEPTH, minv=0, type="pos_numeric", equal=FALSE) || length(filter$MAXDEPTH) != 1){
    warning("Maximum SNP depth filter is invalid. Setting to 500 by defualt.")
    filter$MAXDEPTH <- 500
  }
  if(!is.logical(inferSNPs) || !is.vector(inferSNPs) || length(inferSNPs) != 1)
    stop("Argument for `inferSNPs` not a logical value.")
  
  ## initalize the UR object
  FSobj <- FS$new(RAobj)
  
  ## Define variables that will be used.
  cat("-------------\n")
  cat("Processing Data.\n\n")
  
  ## Extract the private variables we want
  indID <- FSobj$.__enclos_env__$private$indID
  nSnps <- FSobj$.__enclos_env__$private$nSnps
  
  ## check if there are not specific SNPs to keep
  if(is.null(keepSNPs))
    keepSNPs = rep(FALSE, nSnps)
  else if(GUSbase::checkVector(keepSNPs, type="pos_integer", minv=1, maxv=nSnps))
    stop(paste0("Argument for `keepSNPs` is not a vecotr of integers between 0 and the number of SNPs "),nSnps)
  else{
    tmplog = rep(FALSE, nSnps)
    tmplog[keepSNPs] = TRUE
    keepSNPs = tmplog
    rm(tmplog)
  }
  
  ## sort out the pedigree
  ped <- utils::read.csv(pedfile, stringsAsFactors=F)
  famInfo = list()
  ## work out how many families there are
  parents <- unique(ped[c("Mother","Father")])
  # check that parents are in the file
  parents <- parents[which(apply(parents,1, function(x) all(x %in% ped$IndividualID))),]
  ## find the families to be used
  if(is.null(family)){
    family <- unique(ped$Family)[!is.na(unique(ped$Family))]
    family <- family[-which(nchar(family)==0)]
  }
  else{
    if(any(!(family %in% unique(ped$Family)[!is.na(unique(ped$Family))])))
      stop("Family missing from the pedigree file. Please check the family ID suppied or the pedigree file.")
  }
  if(length(family)>1)
    stop("Multiple families are yet to be implemented")
  ## Create each family
  for(fam in family){
    progIndx <- which(ped$Family == fam)
    famParents <- unique( ped[progIndx,c("Mother","Father")] )
    if(nrow(famParents) > 1)
      stop(paste0("Individuals with the same family name (Family ",fam,") have different different parents. Please check the pedigree file"))
    if(any(is.na(famParents)))
      stop(paste0("Family ",fam," has missing parents. Please check the pedigree file"))
    father <- as.integer(famParents["Father"])
    mother <- as.integer(famParents["Mother"])
    ## Create the pedigree
    famInfo[[as.character(fam)]]$progeny <- ped$SampleID[progIndx]
    famInfo[[as.character(fam)]]$parents <- list(Father=ped$SampleID[which(ped$IndividualID==father)],
                                                 Mother=ped$SampleID[which(ped$IndividualID==mother)])
    ## Check to see if there are any grandparents
    # grandparents = unique(ped[which(ped$IndividualID==father),c("Mother","Father")])
    # if(nrow(grandparents) == 1){
    #   famInfo[[as.character(fam)]]$grandparents$paternalGrandFather <- ped$SampleID[which(ped$IndividualID == grandparents[,"Father"])]
    #   famInfo[[as.character(fam)]]$grandparents$paternalGrandMother <- ped$SampleID[which(ped$IndividualID == grandparents[,"Mother"])]
    # }
    # else if(nrow(grandparents) > 1)
    #   stop("Father of family ",fam," has mutiple parents. Please check the pedigree file.")
    # grandparents = unique(ped[which(ped$IndividualID==mother),c("Mother","Father")])
    # if(nrow(grandparents) == 1){
    #   famInfo[[as.character(fam)]]$grandparents$maternalGrandFather <- ped$SampleID[which(ped$IndividualID == grandparents[,"Father"])]
    #   famInfo[[as.character(fam)]]$grandparents$maternalGrandMother <- ped$SampleID[which(ped$IndividualID == grandparents[,"Mother"])]
    # }   
    # else if(nrow(grandparents) > 1)
    #   stop("Mother of family ",fam," has mutiple parents. Please check the pedigree file.")
  }
  
  noFam <- length(famInfo)
  config_all <- config_infer_all <- nInd_all <- indx <- indID_all <- vector(mode = "list", length = noFam)
  genon_all <- ref_all <- alt_all <- vector(mode="list", length=noFam)
  
  ## extract the data and format correct for each family.
  for(fam in 1:noFam){
    cat("Processing Family ",names(famInfo)[fam],".\n\n",sep="")
    mum <- famInfo[[fam]]$parents$Mother
    dad <- famInfo[[fam]]$parents$Father
    # patgrandmum <- famInfo[[fam]]$grandparents$paternalGrandMother
    # patgranddad <- famInfo[[fam]]$grandparents$paternalGrandFather
    # matgrandmum <- famInfo[[fam]]$grandparents$maternalGrandMother
    # matgranddad <- famInfo[[fam]]$grandparents$maternalGrandFather
    ## index the parents
    mumIndx <- which(indID %in% mum)
    if(length(mumIndx) == 0)
      stop(paste0("Mother ID not found family ",fam,"."))
    dadIndx <- which(indID %in% dad)
    if(length(dadIndx) == 0)
      stop(paste0("Father ID not found family ",fam,"."))
    ## index the grandparents
    # patgrandparents <- matgrandparents <- FALSE
    # patgrandmumIndx <- which(indID %in% patgrandmum)
    # patgranddadIndx <- which(indID %in% patgranddad)
    # matgrandmumIndx <- which(indID %in% matgrandmum)
    # matgranddadIndx <- which(indID %in% matgranddad)
    # if((!is.null(patgrandmumIndx) & length(patgrandmumIndx)>0) && (!is.null(patgranddadIndx) & length(patgranddadIndx)>0))
    #   patgrandparents <- TRUE
    # if((!is.null(matgrandmumIndx) & length(matgrandmumIndx)>0) && (!is.null(matgranddadIndx) & length(matgranddadIndx)>0))
    #   matgrandparents <- TRUE
    ## index the progeny
    progIndx <- which(indID %in% famInfo[[fam]]$progeny)
    nInd <- length(progIndx)
    indID_all[[fam]] <- indID[progIndx]
    ## Subset the genon and depth matrices
    genon     <- FSobj$.__enclos_env__$private$genon[progIndx,]
    ref <- FSobj$.__enclos_env__$private$ref[progIndx,]
    alt <- FSobj$.__enclos_env__$private$alt[progIndx,]
    
    ## Determine the segregation types of the loci
    ref_mum = FSobj$.__enclos_env__$private$ref[mumIndx,,drop=FALSE]
    depth_mum = ref_mum + FSobj$.__enclos_env__$private$alt[mumIndx,,drop=FALSE]
    ref_mum = colSums(ref_mum)
    depth_mum = colSums(depth_mum)
    ratio_mum = ref_mum/depth_mum
    genon_mum = (ratio_mum > (1-0.05)) + (ratio_mum > 0.05)
    genon_mum[(depth_mum <= filter$DEPTH) & (genon_mum != 1)] = NA
    
    ref_dad = FSobj$.__enclos_env__$private$ref[dadIndx,,drop=FALSE]
    depth_dad = ref_dad + FSobj$.__enclos_env__$private$alt[dadIndx,,drop=FALSE]
    ref_dad = colSums(ref_dad)
    depth_dad = colSums(depth_dad)
    ratio_dad = ref_dad/depth_dad
    genon_dad = (ratio_dad > (1-0.05)) + (ratio_dad > 0.05)
    genon_dad[(depth_dad <= filter$DEPTH) & (genon_dad != 1)] = NA
    
    # if(patgrandparents){
    #   genon_patgrandmum <- matrix(FSobj$.__enclos_env__$private$genon[patgrandmumIndx,], nrow=length(patgrandmumIndx), ncol=nSnps) 
    #   depth_patgrandmum <- matrix(FSobj$.__enclos_env__$private$ref[patgrandmumIndx,] +
    #                                 FSobj$.__enclos_env__$private$alt[patgrandmumIndx,], nrow=length(patgrandmumIndx), ncol=nSnps)
    #   genon_patgranddad <- matrix(FSobj$.__enclos_env__$private$genon[patgranddadIndx,], nrow=length(patgranddadIndx), ncol=nSnps) 
    #   depth_patgranddad <- matrix(FSobj$.__enclos_env__$private$ref[patgranddadIndx,] +
    #                                 FSobj$.__enclos_env__$private$alt[patgranddadIndx,], nrow=length(patgranddadIndx), ncol=nSnps)
    # }
    # if(matgrandparents){
    #   genon_matgrandmum <- matrix(FSobj$.__enclos_env__$private$genon[matgrandmumIndx,], nrow=length(matgrandmumIndx), ncol=nSnps) 
    #   depth_matgrandmum <- matrix(FSobj$.__enclos_env__$private$ref[matgrandmumIndx,] +
    #                                 FSobj$.__enclos_env__$private$alt[matgrandmumIndx,], nrow=length(matgrandmumIndx), ncol=nSnps)
    #   genon_matgranddad <- matrix(FSobj$.__enclos_env__$private$genon[matgranddadIndx,], nrow=length(matgranddadIndx), ncol=nSnps) 
    #   depth_matgranddad <- matrix(FSobj$.__enclos_env__$private$ref[matgranddadIndx,] +
    #                                 FSobj$.__enclos_env__$private$alt[matgranddadIndx,], nrow=length(matgranddadIndx), ncol=nSnps)
    # }
    
    ## Run the filtering of the progeny SNPs
    MAF <- colMeans(genon, na.rm=T)/2
    MAF[which(is.na(MAF))] = 0.5
    MAF <- pmin(MAF,1-MAF)
    # Check in case SNPs with no genotypes.
    if(any(is.na(MAF))) MAF[is.na(MAF)] = 1
    miss <- apply(genon,2, function(x) sum(is.na(x))/length(x))
    maxdepth <- colMeans(ref + alt)
    indx_temp <- rep(TRUE, nSnps)
    indx_temp[which((MAF < filter$MAF | (miss > filter$MISS) | (miss == 1) | maxdepth > filter$MAXDEPTH) & !keepSNPs)] <- FALSE
    
    ## create file for showing which SNPs fail the filtering threshold.
    SNPfilt = data.frame(CHROM=FSobj$.__enclos_env__$private$chrom, 
                         POS=FSobj$.__enclos_env__$private$pos,
                         SNP_Name=FSobj$.__enclos_env__$private$SNP_Names)
    SNPfilt$MAF = MAF < filter$MAF
    SNPfilt$MISS = (miss > filter$MISS) | (miss == 1)
    SNPfilt$MAXDEPTH = maxdepth > filter$MAXDEPTH

    ## fix max depth to be at 500 total reads
    ## Run into numeric issues otherwise
    high_depth <- which((ref+alt) > 500)
    ref_temp <- ref[high_depth]
    alt_temp <- alt[high_depth]
    ref[high_depth] <- as.integer(ref_temp/(ref_temp+alt_temp)*500)
    alt[high_depth] <- as.integer(alt_temp/(ref_temp+alt_temp)*500)
    
    ## Determine segregation type
    parHap_pat = c("BB","AB","AA")[genon_dad + 1]
    parHap_pat[!indx_temp] = NA
    
    parHap_mat = c("BB","AB","AA")[genon_mum + 1]
    parHap_mat[!indx_temp] = NA
    
    config <- rep(NA,nSnps)
    config[which(parHap_pat == "AB" & parHap_mat == "AB")] <- 1
    config[which(parHap_pat == "AB" & parHap_mat == "AA")] <- 2
    config[which(parHap_pat == "AB" & parHap_mat == "BB")] <- 3
    config[which(parHap_pat == "AA" & parHap_mat == "AB")] <- 4
    config[which(parHap_pat == "BB" & parHap_mat == "AB")] <- 5
    
    SNPfilt$DEPTH = is.na(config)
    SNPfilt$DEPTH[SNPfilt$MAF | SNPfilt$MISS | SNPfilt$MAXDEPTH] = NA
    SNPfilt$PVALUE = FALSE
    SNPfilt$PVALUE[which(SNPfilt$DEPTH | is.na(SNPfilt$DEPTH))] = NA
    
    #### Segregation test to determine if the SNPs have been miss-classified
    seg_Dis <- sapply(1:nSnps,function(x){
      if(is.na(SNPfilt$PVALUE[x]))
        return(NA)
      else{
        d = ref[,x] + alt[,x]
        g = genon[,x]
        K = sum(1/2^(d[which(d != 0)])*0.5)/sum(d != 0)
        nAA = sum(g==2, na.rm=T)
        nAB = sum(g==1, na.rm=T)
        nBB = sum(g==0, na.rm=T)
        ## check that there are sufficient data to perform the chisq test
        if(sum(nAA+nAB+nBB)/length(g) < max(1-filter$MISS,1e-6))
          return(NA)
        else if(config[x] == 1){
          exp_prob <- c(0.25 + K,0.5 - 2*K, 0.25 + K)
          ctest <- suppressWarnings(stats::chisq.test(c(nBB,nAB,nAA), p = exp_prob))
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
        else if(config[x] %in% c(2,4)){
          exp_prob <- c(K, 0.5 - 2*K, 0.5 + K)
          ctest <- suppressWarnings(stats::chisq.test(c(nBB,nAB,nAA), p = exp_prob))
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
        else if(config[x] %in% c(3,5)){
          exp_prob <- c(0.5 + K, 0.5 - 2*K, K)
          ctest <- suppressWarnings(stats::chisq.test(c(nBB,nAB,nAA), p = exp_prob))
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
      }
    },simplify = T)
    config[which(seg_Dis & !keepSNPs)] <- NA
    
    SNPfilt$PVALUE[!is.na(SNPfilt$PVALUE)] = is.na(config[!is.na(SNPfilt$PVALUE)])
    
    ## Infer geotypes for over SNPs that have passed the MAF and MISS thresholds
    #propHeter <- sapply(1:nSnps, function(x) sum(genon[,x] == 1,na.rm=T)/sum(!is.na(genon[,x])))
    config_infer <- rep(NA, nSnps)
    if(inferSNPs){
      toInfer <- which(indx_temp & is.na(config))
      
      seg_Infer <- sapply(toInfer, function(x){
        d = ref[,x] + alt[,x]
        g = genon[,x]
        K = sum(1/2^(d[which(d != 0)])*0.5)/sum(d != 0)
        nAA = sum(g==2, na.rm=T)
        nAB = sum(g==1, na.rm=T)
        nBB = sum(g==0, na.rm=T)
        ## check that there are sufficient data to perform the chisq test
        if(sum(nAA+nAB+nBB)/length(g) < max(1-filter$MISS,1e-6))
          return(NA)
        ## compute chiseq test for both loci types
        exp_prob_BI <- c(0.25 + K,0.5 - 2*K, 0.25 + K)
        exp_prob_SI <- c(K, 0.5 - 2*K, 0.5 + K)
        ctest_BI <- suppressWarnings(stats::chisq.test(c(nBB,nAB,nAA), p = exp_prob_BI))
        ctest_SI_1 <- suppressWarnings(stats::chisq.test(c(nBB,nAB,nAA), p = exp_prob_SI))
        ctest_SI_2 <- suppressWarnings(stats::chisq.test(c(nBB,nAB,nAA), p = rev(exp_prob_SI)))
        ## Check if chisq-test returns NA. If so, don't infer segregation type
        if(any(is.na(ctest_BI$p.value), is.na(ctest_SI_1$p.value), is.na(ctest_SI_2$p.value)))
          return(NA)
        ## do tests to see if we can infer type
        if( ctest_BI$p.value > filter$PVALUE & ctest_SI_1$p.value < filter$PVALUE & ctest_SI_2$p.value < filter$PVALUE )
          return(1)
        else if ( ctest_BI$p.value < filter$PVALUE & ctest_SI_1$p.value > filter$PVALUE & ctest_SI_2$p.value < filter$PVALUE )
          return(4)
        else if ( ctest_BI$p.value < filter$PVALUE & ctest_SI_1$p.value < filter$PVALUE & ctest_SI_2$p.value > filter$PVALUE )
          return(5)
        else
          return(NA)
      },simplify = T)
      config_infer[toInfer] <- seg_Infer
      SNPfilt$PVALUE_infer = is.na(config_infer)
      SNPfilt$PVALUE_infer[which(!SNPfilt$MAF & !SNPfilt$MISS & !SNPfilt$MAXDEPTH & !SNPfilt$DEPTH & !SNPfilt$PVALUE)] = NA
      SNPfilt$PVALUE_infer[which(SNPfilt$MAF | SNPfilt$MISS | SNPfilt$MAXDEPTH )] = NA
    }
    
    chrom <- FSobj$.__enclos_env__$private$chrom
    chrom[(!indx_temp) | (is.na(config) & is.na(config_infer))] <- NA
    pos <- FSobj$.__enclos_env__$private$pos
    pos[(!indx_temp) | (is.na(config) & is.na(config_infer))] <- NA
    ## Extract one SNP from each read.
    set.seed(36475)
    if(filter$BIN > 0){
      oneSNP <- rep(FALSE,nSnps)
      oneSNP[unlist(sapply(unique(chrom), function(x){
        ind <- which(chrom == x)
        if(length(ind) > 0){
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
        } else return(NULL)
      },USE.NAMES = F ))] <- TRUE
      oneSNP[which(keepSNPs)] = TRUE
    } else oneSNP <- rep(TRUE,nSnps)
    
    SNPfilt$BIN = !oneSNP
    if(inferSNPs) SNPfilt$BIN[(is.na(SNPfilt$PVALUE) | SNPfilt$PVALUE) & (SNPfilt$PVALUE_infer | is.na(SNPfilt$PVALUE_infer))] = NA
    else  SNPfilt$BIN[is.na(SNPfilt$PVALUE) | SNPfilt$PVALUE] = NA
    
    indx_temp[which(indx_temp & (!oneSNP | (is.na(config) & is.na(config_infer))))] <- FALSE
    indx[[fam]] <- indx_temp
    config[!indx_temp] <- config_infer[!indx_temp] <- NA
    
    ## Output summary of SNPs filtering:
    cat("Filtering criteria for removing SNPs:\n")
    cat("Minor allele frequency (MAF) < ", filter$MAF,"\n", sep="")
    cat("Percentage of missing genotypes > ", filter$MISS*100,"%\n", sep="")
    cat("Distance for binning SNPs <= ", filter$BIN," bp\n", sep="")
    cat("Read depth associated with at least one parental genotype <= ", filter$DEPTH,"\n", sep="")
    cat("P-value for segregation test < ", filter$PVALUE,"\n", sep="")
    cat("Average SNP depth is > ", filter$MAXDEPTH,"\n\n", sep="")
  
    cat("Filtering SNPs: Number of SNPs failing filtering criteria:\n")
    cat("  MAF:                   ", sum(SNPfilt$MAF & !SNPfilt$MISS & !SNPfilt$MAXDEPTH),"\n", sep="")
    cat("  MISS:                  ", sum(!SNPfilt$MAF & SNPfilt$MISS & !SNPfilt$MAXDEPTH),"\n", sep="")
    cat("  MAXDEPTH:              ", sum(!SNPfilt$MAF & !SNPfilt$MISS & SNPfilt$MAXDEPTH),"\n", sep="")
    cat("  MAF & MISS:            ", sum(SNPfilt$MAF & SNPfilt$MISS & !SNPfilt$MAXDEPTH),"\n", sep="")
    cat("  MAF & MAXDEPTH:        ", sum(SNPfilt$MAF & !SNPfilt$MISS & SNPfilt$MAXDEPTH),"\n", sep="")
    cat("  MISS & MAXDEPTH:       ", sum(!SNPfilt$MAF & SNPfilt$MISS & SNPfilt$MAXDEPTH),"\n", sep="")
    cat("  MAF & MISS & MAXDEPTH: ", sum(SNPfilt$MAF & SNPfilt$MISS & SNPfilt$MAXDEPTH),"\n", sep="")
    cat("  Total discarded:       ",sum(SNPfilt$MAF | SNPfilt$MISS | SNPfilt$MAXDEPTH),"\n",sep="")
    
    cat("Checking remaining SNPs: Number of additional SNPs removed based on:\n")
    cat("  DEPTH:           ",sum(SNPfilt$DEPTH, na.rm=T),"\n",sep="")
    cat("  PVALUE:          ",sum(SNPfilt$PVALUE, na.rm=T),"\n",sep="")
    cat("  Total discarded: ",sum(SNPfilt$PVALUE, na.rm=T) + sum(SNPfilt$DEPTH, na.rm=T),"\n",sep="")
    if(inferSNPs) cat("\nNumber of SNPs with inferred segregation ratio: ",sum(!SNPfilt$PVALUE_infer, na.rm=T),"\n",sep="")
    
    if(filter$BIN > 0) cat("\nA further ",sum(SNPfilt$BIN, na.rm=TRUE)," SNPs were removed based on the BIN filter\n", sep="")
    
    write.csv(SNPfilt, file = "SNP_filtering_GUSMap.csv", row.names=FALSE, quote=FALSE)
    
    ## Determine the segregation groups
    config_all[[fam]] <- config
    config_infer_all[[fam]] <- config_infer
    nInd_all[[fam]] <- nInd
    if(noFam > 1){
      ref[,!indx_temp] <- as.integer(0)
      alt[,!indx_temp] <- as.integer(0)
      genon[,!indx_temp] <- NA
    }
    genon_all[[fam]] <- genon
    ref_all[[fam]] <- ref
    alt_all[[fam]] <- alt
  }
  
  if(noFam == 1){
    fam = 1
    
    ## Find all the SNPs to keep and subset the global variables
    indx_all <- indx[[fam]]
    
    #indx_all <- do.call("rbind",indx)
    #indx_all <- apply(indx_all, 2, any)
    
    genon_all[[fam]]     <- genon[,indx_all]
    ref_all[[fam]] <- ref[,indx_all]
    alt_all[[fam]] <- alt[,indx_all]
    chrom_all            <- FSobj$.__enclos_env__$private$chrom[indx_all]
    pos_all              <- FSobj$.__enclos_env__$private$pos[indx_all]
    SNP_Names            <- FSobj$.__enclos_env__$private$SNP_Names[indx_all]
    
    #group <- group_infer <- vector(mode="list", length=noFam)
    group <- group_infer <- list()
    group$BI <- which(config_all[[fam]][indx_all] == 1)
    group$PI <- which(config_all[[fam]][indx_all] %in% c(2,3))
    group$MI <- which(config_all[[fam]][indx_all] %in% c(4,5))
    
    group_infer$BI <- which(config_infer_all[[fam]][indx_all] == 1) 
    group_infer$SI <- which(config_infer_all[[fam]][indx_all] %in% c(4,5))
    
    config_all[[fam]] <- config_all[[fam]][indx_all]
    config_infer_all[[fam]] <- config_infer_all[[fam]][indx_all]
    
    cat("-------------\n")
    cat("Summary:\n\n")
    cat("Number of SNPs remaining after filtering:",sum(indx_all),"\n")
    if(inferSNPs)
      cat("Number of SNPs with correct segregation type:", sum(!is.na(config_all[[fam]])) ,"\n")
    cat("Both-informative (BI):", length(group$BI),"\n")
    cat("Maternal-informative (MI):", length(group$MI),"\n")
    cat("Paternal-informative (PI):", length(group$PI),"\n")
    if(inferSNPs){
      cat("Number of SNPs with inferred segregation type:", sum(!is.na(config_infer_all[[fam]])),"\n")
      cat("Both-informative (BI):", length(group_infer$BI),"\n")
      cat("Maternal/Paternal-informative (MI or PI):", length(group_infer$SI),"\n")
    }
    cat("Number of progeny:", nInd_all[[fam]],"\n")
    
    ## Create the summary info:
    temp <- ref_all[[1]] + alt_all[[1]]
    summaryInfo <- list()
    if(!inferSNPs){
      summaryInfo$data <- paste0(c(
        "Single Family Linkage analysis:\n\n",
        "Data Summary:\n",
        "Data file:\t", RAobj$.__enclos_env__$private$infilename,"\n",
        "Mean Depth:\t", round(mean(temp),4),"\n",
        "Mean Call Rate:\t",round(sum(temp!=0)/length(temp),4),"\n",
        "Number of ...\n",
        "  Progeny:\t",unlist(nInd),"\n",
        "  MI SNPs:\t",length(group$MI),"\n",
        "  PI SNPs:\t",length(group$PI),"\n",
        "  BI SNPs:\t",length(group$BI),"\n",
        "  Total SNPs:\t",length(unlist(group)),"\n\n"))
    } else{
      summaryInfo$data <- paste0(c(
        "Single Family Linkage analysis:\n\n",
        "Data Summary:\n",
        "Data file:\t", RAobj$.__enclos_env__$private$infilename,"\n",
        "Mean Depth:\t", round(mean(temp),4),"\n",
        "Mean Call Rate:\t",round(sum(temp!=0)/length(temp),4),"\n",
        "Number of ...\n",
        "  Progeny:\t",unlist(nInd),"\n",
        "  SNPs with known segregation\n",
        "    MI SNPs:\t",length(group$MI),"\n",
        "    PI SNPs:\t",length(group$PI),"\n",
        "    BI SNPs:\t",length(group$BI),"\n",
        "    Total SNPs:",length(unlist(group)),"\n",
        "  SNPs with inferred segregation\n",
        "    SI SNPs:\t",length(group_infer$SI),"\n",
        "    BI SNPs:\t",length(group_infer$BI),"\n",
        "    Total SNPs:",length(unlist(group_infer)),"\n\n"))
    }
  }
  else{
    group <- group.temp <- list()
    ## Work out the SNP groupings
    # BI
    tabInf_BI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(1) & indx[[y]]))))
    group.temp$BI <- as.numeric(names(tabInf_BI)[which(tabInf_BI >= MNIF)])
    # MI 
    tabInf_MI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(4,5) & indx[[y]]))))
    group.temp$MI <- setdiff(as.numeric(names(tabInf_MI)[which(tabInf_MI >= MNIF)]), group.temp$BI)
    # PI 
    tabInf_PI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(2,3) & indx[[y]]))))
    group.temp$PI <- setdiff(as.numeric(names(tabInf_PI)[which(tabInf_PI >= MNIF)]), group.temp$BI)
    
    ## Work out which SNPs to keep
    indx_all <- logical(nSnps)
    indx_all[unlist(group.temp)] <- TRUE
    
    add.bi <- group.temp$PI[which(group.temp$PI %in% group.temp$MI)]
    group.temp$BI <- sort(c(group.temp$BI,add.bi))
    group.temp$MI <- setdiff(group.temp$MI,group.temp$BI)
    group.temp$PI <- setdiff(group.temp$PI,group.temp$BI)
    
    group$BI <- which(sort(unique(unlist(group.temp))) %in% group.temp$BI)
    group$MI <- which(sort(unique(unlist(group.temp))) %in% group.temp$MI)
    group$PI <- which(sort(unique(unlist(group.temp))) %in% group.temp$PI)
    
    ## Subset the data
    for(fam in 1:noFam){
      genon_all[[fam]]     <- genon_all[[fam]][,indx_all]
      ref_all[[fam]] <- ref_all[[fam]][,indx_all]
      alt_all[[fam]] <- alt_all[[fam]][,indx_all]
      config_all[[fam]] <- config_all[[fam]][indx_all]
    }
    chrom_all     <- FSobj$.__enclos_env__$private$chrom[indx_all]
    pos_all       <- FSobj$.__enclos_env__$private$pos[indx_all]
    SNP_Names     <- FSobj$.__enclos_env__$private$SNP_Names[indx_all]
    
    cat("-------------\n")
    cat("Summary:\n\n",sep="")
    cat("Number of SNPs remaining after filtering:",sum(indx_all),"\n")
    cat("Both-informative (BI):", length(group$BI),"\n")
    cat("Maternal-informative (MI):", length(group$MI),"\n")
    cat("Paternal-informative (PI):", length(group$PI),"\n")
    cat("Number of progeny in Family...\n")
    for(fam in 1:noFam)
      cat(names(famInfo)[fam],":", nInd_all[[fam]],"\n")
    cat("\n")
    
    ## Create the summary info:
    temp <- do.call("rbind", ref_all) + do.call("list",alt_all)
    summaryInfo <- list()
    summaryInfo$data <- paste0(c(
      "Multiple Family Linkage analysis:\n\n",
      "Data Summary:\n",
      "Data file:\t", RAobj$.__enclos_env__$private$infilename,"\n",
      "Mean Depth:\t", round(mean(temp),4),"\n",
      "Mean Call Rate:\t",round(sum(temp!=0)/length(temp),4),"\n",
      "Number of ...\n",
      "  Progeny:\t",sum(unlist(nInd)),"\n",
      "  MI SNPs:\t",length(group$MI),"\n",
      "  PI SNPs:\t",length(group$PI),"\n",
      "  BI SNPs:\t",length(group$BI),"\n",
      "  Total SNPs:\t",length(unlist(group)),"\n\n"))
    
    group_infer <- NULL
    config_infer_all <- NULL
  }
  
  ## Update the R6 object and return it
  FSobj$.__enclos_env__$private$updatePrivate(list(
    genon = genon_all, ref = ref_all, alt = alt_all, chrom = chrom_all, pos = pos_all,
    group = group, group_infer = group_infer, config_orig = config_all, config_infer_orig = config_infer_all,
    config = config_all, config_infer = config_infer_all, nInd = nInd_all, nSnps = sum(indx_all), noFam = noFam, indID = indID_all, SNP_Names = SNP_Names,
    masked=rep(FALSE,sum(indx_all)), famInfo=famInfo, summaryInfo=summaryInfo)
  )
  
  return(FSobj)
}
