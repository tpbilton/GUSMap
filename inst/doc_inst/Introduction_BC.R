## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set( collapse = TRUE, comment = "#>", eval=TRUE)
library(GUSMap)

## ----simDS, eval = T----------------------------------------------------------
vcffile <- simDS() # extract filename
vcffile            # filename stored in object vcffile

## ----loadData, eval=T---------------------------------------------------------
# convert VCF file to an RA file
rafile <- VCFtoRA(infilename = vcffile$vcf, direct = "./", makePed=TRUE)

## ----rafile, eval=T-----------------------------------------------------------
rafile # file path of RA file

## ----readRA, eval=T-----------------------------------------------------------
RAdata <- readRA(rafile = rafile, sampthres = 0.01, excsamp = NULL)

## ----class_RAdata, eval=T-----------------------------------------------------
class(RAdata)

## ----RAdata, eval=T-----------------------------------------------------------
RAdata

## ----makeBC, eval=T-----------------------------------------------------------
mySpecies <- makeBC(RAobj = RAdata, pedfile = vcffile$ped, inferSNPs = TRUE,
                      filter = list(MAF = 0.05, MISS = 0.5, BIN = 100,
                        DEPTH = 5, PVALUE = 0.01, MAXDEPTH=500))

## ----help_makeBC, eval=F------------------------------------------------------
#  ?`$makeBC`

## ----2pt_rf-------------------------------------------------------------------
mySpecies$rf_2pt(nClust = 3)

## ----save, eval=FALSE---------------------------------------------------------
#  save(mySpecies, file="BCobject.Rdata")

## ----plotChr------------------------------------------------------------------
# plot 2-point rf matrix for BI and MI SNPs
mySpecies$plotChr(mat = "rf", parent = "maternal", lmai=0.5)

## ----plotChr2-----------------------------------------------------------------
# plot 2-point LOD matrix for BI, PI and MI SNPs
mySpecies$plotChr(mat = "LOD", parent = "both", lmai=0.5)

## ----createLG-----------------------------------------------------------------
mySpecies$createLG(parent = "both", LODthres = 5, nComp = 10, reset = FALSE)

## ----help_createLG, eval=F----------------------------------------------------
#  ?`$createLG`

## ----plotLG_mat---------------------------------------------------------------
mySpecies$plotLG(parent = "both", interactive = FALSE, LG = NULL)

## ----createLG_pat-------------------------------------------------------------
mySpecies$createLG(parent = "paternal", LODthres = 10, nComp = 10)

## ----plotLG_pat2--------------------------------------------------------------
mySpecies$plotLG(parent = "both", interactive = F)

## ----removeLG-----------------------------------------------------------------
mySpecies$removeLG(LG = 6)
mySpecies$print(what = "LG-pts")

## ----removeSNP, eval=FALSE----------------------------------------------------
#  mySpecies$removeSNP(snps = c(526,527))

## ----mergeLG------------------------------------------------------------------
mySpecies$mergeLG(LG = c(2,7), mergeTo="paternal")
mySpecies$print(what = "LG-pts")

## ----addSNPs------------------------------------------------------------------
mySpecies$addSNPs(LODthres = 10, nComp = 10)

## ----addBIsnps----------------------------------------------------------------
mySpecies$addBIsnps(LODthres = 10, nComp = 10)

## ----orderLG, , fig.width=8, fig.height=6-------------------------------------
mySpecies$orderLG(mapfun = "haldane", weight = "LOD2", ndim = 30) 

## ----help_orderLG, eval=F-----------------------------------------------------
#  ?`$orderLG`

## ----computeMap---------------------------------------------------------------
mySpecies$computeMap(chrom = NULL, nThreads = 2, rfthres = 0.1)

## ----plotLG-------------------------------------------------------------------
mySpecies$plotLM(LG = NULL, fun = "haldane", col = c(rep("red",3),rep("blue",3)))

## ----plotSyn, fig.width=7, fig.height=7---------------------------------------
mySpecies$plotSyn()

## ----writeLG, eval=F----------------------------------------------------------
#  mySpecies$writeLM(file = "mySpecies")

## ----help_writeLM, eval=F-----------------------------------------------------
#  ?`$writeLM`

