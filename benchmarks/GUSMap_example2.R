# set number of thread from environment variable OMP_NUM_THREADS

library(devtools)
load_all("/projects/nesi00468/chris/gusmap/GUSMap")

num_threads <- max(strtoi(Sys.getenv('OMP_NUM_THREADS')), 1)
cat("Running with num_threads =", num_threads, "\n")

# Code for testing GUSMap

library(GUSMap)
library(onemap)

## Example 2: Large data set:

noChr=2
nSnps=300
noFam=1
set.seed(5721)
config <- list(sapply(1:noChr, function(x) list(sample(c(1,2,4),size=nSnps, prob=c(1,2,2)/5, replace=T)), simplify=T))
simData <- simFS(0.0033,epsilon=0.01,config=config,nInd=200, meanDepth=5, noChr=noChr, seed1=687534, seed2=6772)

simData$rf_2pt(nClust=num_threads)
## plot the results
simData$plotChr(parent="maternal")
simData$plotChr(parent="paternal")

# Solve the likelihood for the hidden Markov model (HMM) and get recombination fraction estimates
simData$rf_est(method="optim")
## Chromosome 1
tmp <- sum(haldane(simData$para$rf_p[[1]]))
cat(tmp, 98.41231, "\n")
tmp <- simData$para$ep[[1]]
cat(tmp, 0.01010545, "\n")
## Chromosome 2
tmp <- sum(haldane(simData$para$rf_p[[2]]))
cat(tmp, 103.2207, "\n")
tmp <- simData$para$ep[[2]]
cat(tmp, 0.009815709, "\n")
