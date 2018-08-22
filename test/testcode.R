

noChr=5
nSnps=50
noFam=1
set.seed(5721)
config <- list(sapply(1:noChr, function(x) list(sample(c(1,2,4),size=nSnps, prob=c(1,2,2)/5, replace=T)), simplify=T))
simData <- simFS(1/nSnps,epsilon=0.01,config=config,nInd=100, meanDepth=5, noChr=noChr, seed1=687534, seed2=6772)

simData$rf_2pt(nClust=3, err=T)
## plot the results
simData$plotChr(parent="maternal")
simData$plotChr(parent="paternal")
simData$plotChr(parent="both")

simData$createLG(parent="both",LODthres = 5)
simData$plotLG(parent="maternal")
simData$plotLG(parent="paternal")
simData$plotLG(parent="both")

simData$addBIsnps(LODthres=5)
simData$plotLG(parent="maternal")
simData$plotLG(parent="paternal")

simData$orderLG(weight = "none", mapfun = "morgan", ndim=30)

simData$rf_est()


simData$setLG()

simData$orderLG(chrom=1, weight = "none", ndim=30)
simData$orderLG(chrom=2, mapfun="morgan", weight = "none", ndim=5)

simData$rf_est()


simData$rf_est(chrom=1, mapped = F)



ind <- c(simData$LG_mat[[1]],simData$LG_pat[[1]])
nSnps <- length(ind)
tt <- expand.grid(1:nSnps, 1:nSnps)
tt <- tt[-which(tt[,1] == tt[,2])]
tempmat <- simData$.__enclos_env__$private$rf[ind,ind]
tt <- cbind(tt,as.vector(tempmat[lower.tri(tempmat)]))

