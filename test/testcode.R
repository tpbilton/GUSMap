
mkfile <- Manuka11()
rafile <- VCFtoRA(mkfile$vcf)
mkdata <- readRA(rafile, gform="reference")
mkfs <- makeFS(mkdata, pedfile = mkfile$ped)


### Simulated data
noChr=5
nSnps=50
noFam=1
set.seed(5721)
config <- list(sapply(1:noChr, function(x) list(sample(c(1,2,4),size=nSnps, prob=c(1,2,2)/5, replace=T)), simplify=T))
simData <- simFS(1/nSnps,epsilon=0.01,config=config,nInd=100, meanDepth=5, noChr=noChr, seed1=687534, seed2=6772)

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
