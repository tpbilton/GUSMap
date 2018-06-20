# set number of thread from environment variable OMP_NUM_THREADS

num_threads <- max(strtoi(Sys.getenv('OMP_NUM_THREADS')), 1)
cat("Running with num_threads =", num_threads, "\n")

# Code for testing GUSMap

library(GUSMap)
library(onemap)

## Example 1: Manuka dataset 

#### Read in data

## get file paths to Manuka data
MKdata <- Manuka11()
## Convert VCF to RA
RAfile <- VCFtoRA(MKdata$vcf)
## read in the RA fill
RAdata <- readRA(RAfile, gform = "reference")

#### Create a mapping population

MK_fs <- createPop(RAdata, pedfile=MKdata$ped)

#### Compute pairwise rf's

MK_fs$rf_2pt(nClust=num_threads)

MK_fs$plotChr(parent = "maternal")
MK_fs$plotChr(parent = "paternal")

#### Mask SNPs

##Want to mask the SNPs know to be bad
badSnps <- rep(FALSE, MK_fs$.__enclos_env__$private$nSnps)
badSnps[c(3:14,16,18:23,25:27,29,33:35,39:42,44,46:49, 
   53,58:60,63,65:67,69:74,76,77,85,104,114:129, 
   152,162,163,173:180,187:190,206:208,210,212,213, 
   216,228,272,285,286,298,304,314,318,320,342,344,345, 
   354,358,370,373,380,389,390,397,408,411,421,427, 
   437,446,469,493,498:500,506,510,516,519,522,531, 
   543,553,556,569,582,585,601,617,620,623,628,635, 
   636:643,647:650,669:680)] <- TRUE 
MK_fs$maskSNP(which(badSnps))
MK_fs$plotChr()
MK_fs$plotChr(parent = "paternal")

#### Compute the linkage maps for SNPs on chromosome 11

depth <- MK_fs$.__enclos_env__$private$ref[[1]] + MK_fs$.__enclos_env__$private$alt[[1]]
rd6 <- ((colMeans(depth) > 6) | badSnps)
MK_fs$maskSNP(which(rd6))
## Compute the rf estimates
## using optim 
MK_fs$rf_est(method="optim")
tmp <- sum(haldane(MK_fs$para$rf_p[[1]])) ## should be 76.09481
cat(tmp, 76.09481, "\n")
tmp <- MK_fs$para$ep[[1]]                 ## should be 0.003167973
cat(tmp, 0.003167973, "\n")
## using EM algorithm
MK_fs$rf_est(method="EM")
tmp <- sum(haldane(MK_fs$para$rf_p[[1]])) ## should be 76.09481
cat(tmp, 76.09481, "\n")
tmp <- MK_fs$para$ep[[1]]                 ## should be 0.003167973
cat(tmp, 0.003167973, "\n")
