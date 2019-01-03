
mkfile <- Manuka11()
rafile <- VCFtoRA(mkfile$vcf)
mkdata <- readRA(rafile, gform="reference")
mkfs <- makeFS(mkdata, pedfile = mkfile$ped)


### Simulated data
noChr=1
nSnps=50
noFam=1
set.seed(5721)
config <- list(sapply(1:noChr, function(x) list(sample(c(1,2,4),size=nSnps, prob=c(1,2,2)/5, replace=T)), simplify=T))
simData <- simFS(0.001,epsilon=0.0001,config=config,nInd=100, meanDepth=25, seed1=687534, seed2=6772)

simData$rf_2pt(nClust=3, err=T)
## plot the results
simData$plotChr(parent="maternal", lmai = 0.5)
simData$plotChr(parent="paternal", lmai = 0.5)
simData$plotChr(parent="both", lmai=0.5)

simData$createLG(parent="both",LODthres = 5)
simData$plotLG(parent="maternal")
simData$plotLG(parent="paternal")
simData$plotLG(parent="both")

simData$addBIsnps(LODthres=5)
simData$plotLG(parent="maternal")
simData$plotLG(parent="paternal")

simData$orderLG(weight = "none", mapfun = "morgan", ndim=30)

simData$rf_est()

simData$plotSyn()


## test viterbi algorithm
simData$computeMap(mapped = FALSE)
res <- simData$extractVar("para")
nSnps <- simData$extractVar("nSnps")[[1]]
nInd <- simData$extractVar("nInd")[[1]][[1]]

tt <- viterbi_fs_err(res$para$rf_p[[1]], rep(res$para$ep[[1]][[1]],nSnps), nInd, nSnps,
                     simData$extractVar("ref")[[1]][[1]], simData$extractVar("alt")[[1]][[1]], res$para$OPGP[[1]])
infer_m <- rbind(tt < 2, apply(tt, 2, function(x) x %in% c(0,2))) + 0
newrf <- sapply(2:nSnps, function(x) sum(infer_m[,x]!=infer_m[,x-1]))/(2*nInd)
## if we want to remove stuff where recombination is not informative
infer_m[1:100,which(res$para$OPGP[[1]] %in% c(9:12))] <- NA
infer_m[101:200,which(res$para$OPGP[[1]] %in% c(5:8))] <- NA
newrf <- sapply(2:nSnps, function(x) sum(infer_m[,x]!=infer_m[,x-1], na.rm = T))

temp <- sapply(2:nSnps, function(x) sapply(1:(2*nInd), function(y) {
  tt <- is.na(infer_m[y,1:(x-1)])
  if(all(tt))
     infer_m[y,x] = FALSE
  else
    infer_m[y,x]!= infer_m[y,tail(which(!tt),1)]
}))
newrf <- apply(temp, 2, function(x) sum(x, na.rm=T)/length(!is.na(x)))

infer_m[1:nInd,OPGP %in% c(9:12)] <- NA
infer_m[1:nInd+nInd], OPGP %in% c(5:8)] <- NA
sum(GUSMap:::mfun(newrf/(2*nInd), centiM=T))




