### R Function for reading in the Manuka data of chromosome 11.
### Author: Timothy Bilton
### Date: 29/08/17
### Edited: 18/10/17


## Wrapper function for reading in the Manuka data set for chromosome 11
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
  indID <- switch(gform, uneak = ghead[2:(nInd + 1)], Tassel = ghead[-(1:2)])
  
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
