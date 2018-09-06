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
#' Make a full-sib family (FS) population
#'
#' Create an FS object from an RA object, perform standard filtering and computes statistics specific to full-sib family populations.
#' 
#' The segregation type of each SNP is infer based on the genotypes of the parents. The parental genotypes are called homozygous for the 
#' reference is there is only reference reads seen, heterozygous if at least one read for the reference and alternate allele are seen
#' and homozygous for the alternate if only reads for the alternate allele is seen. 
#' The \code{BIN} filter is implemented to remove any SNPs
#' with incorrect segregation type. The segregation test procedure for performing the segregation test is described in the supplementary methods 
#' of the publication by \insertCite{bilton2018genetics1;textual}{GUSMap} (Section 4 of File S1).
#' 
#' @param RAobj Object of class RA created via the \code{\link[GUSbase]{readRA}} function.
#' @param pedfile Character string giving the file name (relative to the current directory) of the pedigree file.
#' @param family Vector of character strings giving the families to retain in the FS object. This allows a pedigree file with more than one family to be supplied.
#' @param filter Named list of thresholds for various filtering criteria.
#' See below for details.
#' 
#' @section Filtering:
#' The filtering criteria currently implemented are:
#' \itemize{
#' \item{Minor allele frequency (\code{MAF}): }{SNPs are discarded if their MAF is less than the threshold (default is 0.05)}
#' \item{Proportion of missing data (\code{MISS}): }{SNPs are discarded if the proportion of individuals with no reads (e.g. missing genotype)
#'  is greater than the threshold value (default is 0.5).}
#' \item{Bin size for SNP selection (\code{BIN}): }{SNPs are binned together if the distance between them is less than the threshold value of base pairs (default is 0).
#' One SNP is then randomly selected from each bin and retained for final analysis. This filtering is to ensure that there is only one SNP on each tag.}
#' \item{Parental read depth (\code{DEPTH}): }{SNPs are discarded if the read depth of either parent is less than the threshold value (default is 6). 
#' This filter is to remove SNPs where the parental information is insufficient to infer segregation type accurately.}
#' \item{Segregation test P-value (\code{PVALUE}): }{SNPs are discarded if the p-value from a segregation test is smaller than the threshold (default is 0.01).
#'  This is filter out SNPs where the segregation type has been inferred wrong.}
#' }
#' 
#' @section Pedigree File:
#' The pedigree file must be a csv file contains the five columns:
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
#' @return 
#' An R6 object of class RA
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

### Make a full-sib family population
makeFS <- function(RAobj, pedfile, family=NULL, 
                   filter=list(MAF=0.05, MISS=0.2, BIN=0, DEPTH=5, PVALUE=0.01)){
  #inferSNPs = FALSE, perInfFam=1){
  inferSNPs = FALSE; perInfFam=1 # some variables for multiple familes needed for later
  ## Do some checks
  if(!all(class(RAobj) %in% c("RA","R6")))
    stop("First argument supplied is not of class 'R6' and 'RA'")
  if(is.null(filter$MAF) || filter$MAF<0 || filter$MAF>1 || !is.numeric(filter$MAF)){
    warning("Minor allele filter has not be specified or is invalid. Setting to 5%:")
    filter$MAF <- 0.05
  }
  if(perInfFam <= 0.5 || perInfFam > 1)
    stop("The of the percentage of families which are informative for each SNP must greater than 50% and leass than equal to 100%.")
  if(is.null(filter$MISS) || filter$MISS<0 || filter$MISS>1 || !is.numeric(filter$MISS)){
    warning("Proportion of missing data filter has not be specified or is invalid. Setting to 20%:")
    filter$MISS <- 0.2
  }
  if(is.null(filter$BIN) || filter$BIN<0 || !is.finite(filter$BIN) || !is.numeric(filter$BIN)){
    warning("Minimum distance between adjacent SNPs is not specified or is invalid. Setting to 0:")
    filter$BIN <- 0 
  }
  if(is.null(filter$DEPTH) || filter$DEPTH<0 || is.infinite(filter$DEPTH) || !is.numeric(filter$DEPTH)){
    warning("Minimum depth on the parental genotypes filter has not be specified or is invalid. Setting to a depth of 5")
    filter$DEPTH <- 5
  }
  if(is.null(filter$PVALUE) || filter$PVALUE<0 || is.infinite(filter$PVALUE) || !is.numeric(filter$PVALUE)){
    warning("P-value for segregation test is not specified or invalid. Setting a P-value of 0.01:")
    filter$PVALUE <- 0.05
  }
  
  ## initalize the UR object
  FSobj <- FS$new(RAobj)
  
  ## Define variables that will be used.
  cat("-------------\n")
  cat("Processing Data.\n\n")
  
  cat("Filtering criteria for removing SNPs:\n")
  cat("Minor allele frequency (MAF) < ", filter$MAF,"\n", sep="")
  cat("Percentage of missing genotypes > ", filter$MISS*100,"%\n", sep="")
  cat("Distance for binning SNPs <= ", filter$BIN," bp\n", sep="")
  cat("Read depth associated with at least one parental genotype <= ", filter$DEPTH,"\n", sep="")
  cat("P-value for segregation test < ", filter$PVALUE,"%\n\n", sep="")
  ## Extract the private variables we want
  indID <- FSobj$.__enclos_env__$private$indID
  nSnps <- FSobj$.__enclos_env__$private$nSnps
  
  ## sort out the pedigree
  ped <- read.csv(pedfile, stringsAsFactors=F)
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
    stop("Multiple are yet to be implemented")
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
    grandparents = ped[which(ped$IndividualID==father),c("Mother","Father")]
    if(nrow(grandparents) == 1){
      famInfo[[as.character(fam)]]$grandparents$paternalGrandFather <- ped$SampleID[which(ped$IndividualID == grandparents[,"Father"])]
      famInfo[[as.character(fam)]]$grandparents$paternalGrandMother <- ped$SampleID[which(ped$IndividualID == grandparents[,"Mother"])]
    }
    else if(nrow(grandparents) > 2)
      stop("Father of family ",fam," has mutiple parents. Please check the pedigree file.")
    grandparents = ped[which(ped$IndividualID==mother),c("Mother","Father")]
    if(nrow(grandparents) == 1){
      famInfo[[as.character(fam)]]$grandparents$maternalGrandFather <- ped$SampleID[which(ped$IndividualID == grandparents[,"Father"])]
      famInfo[[as.character(fam)]]$grandparents$maternalGrandMother <- ped$SampleID[which(ped$IndividualID == grandparents[,"Mother"])]
    }   
    else if(nrow(grandparents) > 2)
      stop("Mother of family ",fam," has mutiple parents. Please check the pedigree file.")
  }
  
  noFam <- length(famInfo)
  config_all <- config_infer_all <- nInd_all <- indx <- indID_all <- vector(mode = "list", length = noFam)
  genon_all <- ref_all <- alt_all <- vector(mode="list", length=noFam)
  
  ## extract the data and format correct for each family.
  for(fam in 1:noFam){
    cat("Processing Family ",names(famInfo)[fam],".\n\n",sep="")
    mum <- famInfo[[fam]]$parents$Mother
    dad <- famInfo[[fam]]$parents$Father
    patgrandmum <- famInfo[[fam]]$grandparents$paternalGrandMother
    patgranddad <- famInfo[[fam]]$grandparents$paternalGrandFather
    matgrandmum <- famInfo[[fam]]$grandparents$maternalGrandMother
    matgranddad <- famInfo[[fam]]$grandparents$maternalGrandFather
    ## index the parents
    mumIndx <- which(indID %in% mum)
    if(length(mumIndx) == 0)
      stop(paste0("Mother ID not found family ",fam,"."))
    dadIndx <- which(indID %in% dad)
    if(length(dadIndx) == 0)
      stop(paste0("Father ID not found family ",fam,"."))
    ## index the grandparents
    patgrandparents <- matgrandparents <- FALSE
    patgrandmumIndx <- which(indID %in% patgrandmum)
    patgranddadIndx <- which(indID %in% patgranddad)
    matgrandmumIndx <- which(indID %in% matgrandmum)
    matgranddadIndx <- which(indID %in% matgranddad)
    if(!is.null(patgrandmumIndx) && !is.null(patgranddadIndx))
      patgrandparents <- TRUE
    if(!is.null(matgrandmumIndx) && !is.null(matgranddadIndx))
      matgrandparents <- TRUE
    ## index the progeny
    progIndx <- which(indID %in% famInfo[[fam]]$progeny)
    nInd <- length(progIndx)
    indID_all[[fam]] <- indID[progIndx]
    ## Subset the genon and depth matrices
    genon     <- FSobj$.__enclos_env__$private$genon[progIndx,]
    ref <- FSobj$.__enclos_env__$private$ref[progIndx,]
    alt <- FSobj$.__enclos_env__$private$alt[progIndx,]
    
    ## Determine the segregation types of the loci
    genon_mum <- matrix(FSobj$.__enclos_env__$private$genon[mumIndx,], nrow=length(mumIndx), ncol=nSnps) 
    genon_dad <- matrix(FSobj$.__enclos_env__$private$genon[dadIndx,], nrow=length(mumIndx), ncol=nSnps)
    depth_mum <- matrix(FSobj$.__enclos_env__$private$ref[mumIndx,] +
                          FSobj$.__enclos_env__$private$alt[mumIndx,], nrow=length(mumIndx), ncol=nSnps)
    depth_dad <- matrix(FSobj$.__enclos_env__$private$ref[dadIndx,] +
                          FSobj$.__enclos_env__$private$alt[dadIndx,], nrow=length(mumIndx), ncol=nSnps)
    
    if(patgrandparents){
      genon_patgrandmum <- matrix(FSobj$.__enclos_env__$private$genon[patgrandmumIndx,], nrow=length(patgrandmumIndx), ncol=nSnps) 
      depth_patgrandmum <- matrix(FSobj$.__enclos_env__$private$ref[patgrandmumIndx,] +
                                    FSobj$.__enclos_env__$private$alt[patgrandmumIndx,], nrow=length(patgrandmumIndx), ncol=nSnps)
      genon_patgranddad <- matrix(FSobj$.__enclos_env__$private$genon[patgranddadIndx,], nrow=length(patgranddadIndx), ncol=nSnps) 
      depth_patgranddad <- matrix(FSobj$.__enclos_env__$private$ref[patgranddadIndx,] +
                                    FSobj$.__enclos_env__$private$alt[patgranddadIndx,], nrow=length(patgranddadIndx), ncol=nSnps)
    }
    if(matgrandparents){
      genon_matgrandmum <- matrix(FSobj$.__enclos_env__$private$genon[matgrandmumIndx,], nrow=length(matgrandmumIndx), ncol=nSnps) 
      depth_matgrandmum <- matrix(FSobj$.__enclos_env__$private$ref[matgrandmumIndx,] +
                                    FSobj$.__enclos_env__$private$alt[matgrandmumIndx,], nrow=length(matgrandmumIndx), ncol=nSnps)
      genon_matgranddad <- matrix(FSobj$.__enclos_env__$private$genon[matgranddadIndx,], nrow=length(matgranddadIndx), ncol=nSnps) 
      depth_matgranddad <- matrix(FSobj$.__enclos_env__$private$ref[matgranddadIndx,] +
                                    FSobj$.__enclos_env__$private$alt[matgranddadIndx,], nrow=length(matgranddadIndx), ncol=nSnps)
    }
    
    parHap_pat <- sapply(1:nSnps,function(x){
      x_p = genon_dad[,x]; d_p = depth_dad[,x]
      if(any(x_p==1,na.rm=T))
        return("AB")
      else if(sum(d_p) > filter$DEPTH){
        if(all(x_p==2, na.rm=T))
          return("AA")
        else if(all(x_p==0, na.rm=T))
          return("BB")
        else if(patgrandparents){
          if(sum(depth_patgranddad[,x])>filter$DEPTH && sum(depth_patgrandmum)>filter$DEPTH){
            x_gp = genon_patgranddad[,x]; x_gm = genon_patgrandmum[,x]
            if((x_gp == 2 & x_gm == 0) || (x_gp == 0 & x_gm == 2))
              return("AB")
            else if(x_gp == 2 & x_gm == 2)
              return("AA")
            else if(x_gp == 0 & x_gm == 0)
              return("BB")
          }
        }
        else
          return(NA)
      }
      else
        return(NA)
    })
    
    parHap_mat <- sapply(1:nSnps,function(x){
      x_m = genon_mum[,x]; d_m = depth_mum[,x]
      if(any(x_m==1,na.rm=T))
        return("AB")
      else if(sum(d_m) > filter$DEPTH){
        if(all(x_m==2, na.rm=T))
          return("AA")
        else if(all(x_m==0, na.rm=T))
          return("BB")
        else if(matgrandparents){
          if(sum(depth_matgranddad[,x])>filter$DEPTH && sum(depth_matgrandmum)>filter$DEPTH){
            x_gp = genon_matgranddad[,x]; x_gm = genon_matgrandmum[,x]
            if((x_gp == 2 & x_gm == 0) || (x_gp == 0 & x_gm == 2))
              return("AB")
            else if(x_gp == 2 & x_gm == 2)
              return("AA")
            else if(x_gp == 0 & x_gm == 0)
              return("BB")
          }
        }
        else
          return(NA)
      }
      else
        return(NA)
    })
    
    config <- rep(NA,nSnps)
    config[which(parHap_pat == "AB" & parHap_mat == "AB")] <- 1
    config[which(parHap_pat == "AB" & parHap_mat == "AA")] <- 2
    config[which(parHap_pat == "AB" & parHap_mat == "BB")] <- 3
    config[which(parHap_pat == "AA" & parHap_mat == "AB")] <- 4
    config[which(parHap_pat == "BB" & parHap_mat == "AB")] <- 5
    
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
          ctest <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = exp_prob))
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
        else if(config[x] %in% c(2,4)){
          exp_prob <- c(K, 0.5 - 2*K, 0.5 + K)
          ctest <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = exp_prob))
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
        else if(config[x] %in% c(3,5)){
          exp_prob <- c(0.5 + K, 0.5 - 2*K, K)
          ctest <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = exp_prob))
          return(ifelse(ctest$p.value < filter$PVALUE, TRUE, FALSE))
        }
      }
    },simplify = T)
    config[which(seg_Dis)] <- NA
    
    ## Run the filtering of the progeny SNPs
    MAF <- colMeans(genon, na.rm=T)/2
    MAF <- pmin(MAF,1-MAF)
    miss <- apply(genon,2, function(x) sum(is.na(x))/length(x))
    
    ## Infer geotypes for over SNPs that have passed the MAF and MISS thresholds
    #propHeter <- sapply(1:nSnps, function(x) sum(genon[,x] == 1,na.rm=T)/sum(!is.na(genon[,x])))
    if(inferSNPs){
      toInfer <- (MAF > filter$MAF) & (miss < filter$MISS) & is.na(config)
      
      seg_Infer <- sapply(1:nSnps, function(x){
        if(!toInfer[x])
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
          ## compute chiseq test for both loci types
          exp_prob_BI <- c(0.25 + K,0.5 - 2*K, 0.25 + K)
          exp_prob_SI <- c(K, 0.5 - 2*K, 0.5 + K)
          ctest_BI <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = exp_prob_BI))
          ctest_SI_1 <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = exp_prob_SI))
          ctest_SI_2 <- suppressWarnings(chisq.test(c(nBB,nAB,nAA), p = rev(exp_prob_SI)))
          ## do tests to see if we can infer type
          if( ctest_BI$p.value > filter$PVALUE & ctest_SI_1$p.value < filter$PVALUE & ctest_SI_2$p.value < filter$PVALUE )
            return(1)
          else if ( ctest_BI$p.value < filter$PVALUE & ctest_SI_1$p.value > filter$PVALUE & ctest_SI_2$p.value < filter$PVALUE )
            return(4)
          else if ( ctest_BI$p.value < filter$PVALUE & ctest_SI_1$p.value < filter$PVALUE & ctest_SI_2$p.value > filter$PVALUE )
            return(5)
          else
            return(NA)
        }
      },simplify = T)
    }
    
    chrom <- FSobj$.__enclos_env__$private$chrom
    pos <- FSobj$.__enclos_env__$private$pos
    ## Extract one SNP from each read.
    if(filter$BIN > 0){
      oneSNP <- rep(FALSE,nSnps)
      oneSNP[unlist(sapply(unique(chrom), function(x){
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
    
    if(inferSNPs){
      indx[[fam]] <- (MAF > filter$MAF) & (miss < filter$MISS) & ( !is.na(config) | !is.na(seg_Infer) ) & oneSNP
      config[!indx[[fam]]] <- seg_Infer[!indx[[fam]]] <- NA
    }
    else
      indx[[fam]] <- (MAF > filter$MAF) & (miss < filter$MISS) & (!is.na(config)) & oneSNP
    
    
    ## Determine the segregation groups
    config_all[[fam]] <- config
    if(inferSNPs)
      config_infer_all[[fam]] <- seg_Infer
    
    nInd_all[[fam]] <- nInd
    
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
    if(inferSNPs)
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
  }
  else{
    noInfoFam <- ceiling(perInfFam*noFam)
    group <- group.temp <- list()
    ## Work out the SNP groupings
    # BI
    tabInf_BI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(1) & indx[[y]]))))
    group.temp$BI <- as.numeric(names(tabInf_BI)[which(tabInf_BI >= noInfoFam)])
    # MI 
    tabInf_MI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(4,5) & indx[[y]]))))
    group.temp$MI <- setdiff(as.numeric(names(tabInf_MI)[which(tabInf_MI >= noInfoFam)]), group.temp$BI)
    # PI 
    tabInf_PI <- table(unlist(sapply(1:noFam,function(y) which(config_all[[y]] %in% c(2,3) & indx[[y]]))))
    group.temp$PI <- setdiff(as.numeric(names(tabInf_PI)[which(tabInf_PI >= noInfoFam)]), group.temp$BI)
    
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
    
    group_infer <- NULL
    config_infer_all <- NULL
    summaryInfo  <- NULL
  }
  
  ## Update the R6 object and return it
  FSobj$.__enclos_env__$private$updatePrivate(list(
    genon = genon_all, ref = ref_all, alt = alt_all, chrom = chrom_all, pos = pos_all,
    group = group, group_infer = group_infer, config = config_all, config_infer = config_infer_all,
    nInd = nInd_all, nSnps = sum(indx_all), noFam = noFam, indID = indID_all, SNP_Names = SNP_Names,
    masked=rep(FALSE,sum(indx_all)), famInfo=famInfo, summaryInfo=summaryInfo)
  )
  
  return(FSobj)
}
