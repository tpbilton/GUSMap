### R Function for reading in the Manuka data of chromosome 11.
### Author: Timothy Bilton
### Date: 29/08/17
### Edited: 18/10/17


## Wrapper function for reading in the Manuka data set for chromosome 11


#' Manuka chromosome 11 SNPs
#' 
#' Function for extracting the Manuka data used in the publication by Bilton et
#' al. (2017).
#' 
#' The data consists of 680 SNPs, genotyped using genotyping-by-sequencing
#' methods. The original data is in RA format and is converted into matrices of
#' allele counts for the reference and alternate allele as required by GUSMap.
#' 
#' @return Function outputs a list containing the following elements:
#' \itemize{
#' \item genon Numeric matrix of genotype calls. 2=homozygous for reference,
#' 1=heterozygous, 0=homozygous for alternate and NA missing genotype
#' \item depth_Ref: Numeric matrix of allele counts for the reference allele.
#' \item depth_Alt: Numeric matrix of allele counts for the alternate allele.
#' \item indNames: Names of the individuals. The \emph{i}th entry corresponds
#' to the \emph{i}th row of the \code{genon} and \code{depth} matrices.
#' \item chrom: Chromosome number of the SNPs. The \emph{i}th entry corresponds
#' to the \emph{i}th column of the \code{genon} and \code{depth} matrices.
#' \item pos: Position number of the SNPs. The \emph{i}th entry corresponds to
#' the \emph{i}th column of the \code{genon} and \code{depth} matrices.
#' }
#' @author Timothy P. Bilton
#' @references Bilton, T.P., Schofield, M.R., Black, M.A., Chagne D., Wilcox
#' P., Dodds K.G. (2017). Accounting for errors in low coverage sequencing data
#' when constructing genetic maps using multiparental outcrossed populations.
#' Unpublished manuscript.
#' @examples
#' 
#' ## To convert the data into the required format for GusMap:
#' manuka <- Manuka11()
#' 
#' ## To access the RA directly
#' RAname <- system.file("extdata", "Manuka11.txt.gz", package="GUSMap")
#' RAdata <- read.table(RAname, header=T)
#' 
#' ## Work out the segregation types
#' config <- unlist(apply(manuka$genon[1:4,],2,function(x){
#'   x_p = x[1:2]; x_m = x[3:4]
#'   if(any(x_p==1,na.rm=T) & any(x_m==1,na.rm=T))
#'     return(1)
#'   else if(all(x_m==2,na.rm=T) & any(x_p==1, na.rm=T))
#'     return(2)
#'   else if(all(x_m==0,na.rm=T) & any(x_p==1, na.rm=T))
#'     return(3)
#'   else if(all(x_p==2,na.rm=T) & any(x_m==1,na.rm=T))
#'     return(4)
#'   else if(all(x_p==0,na.rm=T) & any(x_m==1,na.rm=T))
#'     return(5)
#' }))
#' 
#' ## SNPs removed due to mis-ordering
#' badSnps <- rep(FALSE, ncol(manuka$genon))
#' badSnps[c(3:14,16,18:23,25:27,29,33:35,39:42,44,46:49,
#'   53,58:60,63,65:67,69:74,76,77,85,104,114:129,
#'   152,162,163,173:180,187:190,206:208,210,212,213,
#'   216,228,272,285,286,298,304,314,318,320,342,344,345,
#'   354,358,370,373,380,389,390,397,408,411,421,427,
#'   437,446,469,493,498:500,506,510,516,519,522,531,
#'   543,553,556,569,582,585,601,617,620,623,628,635,
#'   636:643,647:650,669:680)] <- TRUE
#' 
#' ## Compute the map for the low coverage set of SNPs 
#' ind_rd6 <- which(colMeans(manuka$depth_Ref[-c(1:4),] + manuka$depth_Alt[-c(1:4),],na.rm=T ) < 6 & !badSnps) 
#' OPGP_rd6 <- infer_OPGP_FS(manuka$depth_Ref[-c(1:4),ind_rd6], manuka$depth_Alt[-c(1:4),ind_rd6],
#'                           config[ind_rd6], epsilon=0.001, reltol=1e-3)
#' 
#' rf_est <- rf_est_FS(init_r=0.01,depth_Ref=list(manuka$depth_Ref[-c(1:4),ind_rd6]),
#'                                 depth_Alt=list(manuka$depth_Alt[-c(1:4),ind_rd6]),
#'                                 OPGP=list(OPGP_rd6), epsilon=0.001)
#' 
#' @export Manuka11
Manuka11 <- function(){
  ## Set the parameters required for KGD
  genofile <- system.file("extdata", "Manuka11.txt.gz", package="GUSMap")
  return(RAtoCount(genofile,"Tassel"))
}


## Function for converting RA data to allele count format
RAtoCount <- function(genofile, gform){
  
  ## separate character between reference and alternate allele count
  gsep <- switch(gform, uneak = "|", Tassel = ",")
  
  ## Process the individuals info
  ghead <- scan(genofile, what = "", nlines = 1, sep = "\t")
  nInd <- length(ghead) - switch(gform, uneak = 6, Tassel = 2)
  indID <- switch(gform, uneak = ghead[2:(nind + 1)], Tassel = ghead[-(1:2)])
  
  ## Read in the data
  # If reference based
  if (gform == "Tassel"){
    genosin <- scan(genofile, skip = 1, sep = "\t", what = c(list(chrom = "", coord = 0), rep(list(""), nInd)))
    chrom <- genosin[[1]]
    pos <- genosin[[2]]
    SNP_Names <- paste(genosin[[1]],genosin[[2]],sep="_")
  }
  # Non-reference based 
  if (gform == "uneak"){
    genosin <- scan(genofile, skip = 1, sep = "\t", what = c(list(chrom = ""), rep(list(""), nind), list(hetc1 = 0, hetc2 = 0, acount1 = 0, acount2 = 0, p = 0)))
    SNP_Names <- genosin[[1]]
  }
  ## Compute the number of SNPs
  nSnps <- length(SNP_Names)
  
  ## Compute the read count matrices and genon matrix 
  depth_Ref <- depth_Alt <- matrix(0, nrow = nInd, ncol = nSnps)
  for (i in 1:nInd){ 
    depths <- strsplit(genosin[[i + switch(gform, uneak = 1, Tassel = 2)]], split = gsep, fixed = TRUE)
    depth_Ref[i, ] <- as.numeric(unlist(lapply(depths,function(z) z[1])))
    depth_Alt[i, ] <- as.numeric(unlist(lapply(depths,function(z) z[2])))
  }
  genon <- (depth_Ref > 0) + (depth_Alt == 0)
  genon[which(depth_Ref == 0 & depth_Alt == 0)] <- NA
  if (gform == "uneak") 
    AFrq <- genosin[[length(genosin)]]
  
  ## return an object of the type we want
  return(
    list(genon=genon, depth_Ref=depth_Ref,depth_Alt=depth_Alt,chrom=chrom,pos=pos,indID=indID,SNP_Names=SNP_Names)
  )
}
