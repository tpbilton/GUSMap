### R Function for reading in the Manuka data of chromosome 11.
### Author: Timothy Bilton
### Date: 29/08/17

Manuka11 <- function(){
  ## Set the parameters required for KGD
  genofile <- system.file("extdata", "Manuka11.txt.gz", package="GUSMap")
  gform <- "Tassel"
  ## Run KGD to obtain the genon, depth and seqID objects
  source("https://raw.github.com/AgResearch/KGD/master/GBS-Chip-Gmatrix.R", local=environment())
  toReturn <- list(genon=genon,depth=depth.orig,indNames=seqID,chrom=chrom,pos=pos)
  ## Remove the objects created unnecessarily.
  rm(list=c("allelecounts","alleles","chrom","depth","depth.orig","fcolo","genon","HWdis","pg","pos","RAcounts","sampdepth","l10LRT","maf","nind","nsnps","p","sampdepth.max","sampdepth.med","samples","seqID","SNP_Names","snpdepth"), envir = globalenv())
  return(toReturn)
}