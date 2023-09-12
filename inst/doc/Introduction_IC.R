## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set( collapse = TRUE, comment = "#>", eval=TRUE)
library(GUSMap)

## ----simDS, eval = T----------------------------------------------------------
vcffile <- simDS() # extract filename
vcffile            # filename stored in object vcffile

## ----loadData, eval=T---------------------------------------------------------
# convert VCF file to an RA file
rafile <- VCFtoRA(infilename = vcffile$vcf, direct = "./")

## ----rafile, eval=T-----------------------------------------------------------
rafile # file path of RA file

## ----readRA, eval=T-----------------------------------------------------------
RAdata <- readRA(rafile = rafile, sampthres = 0.01, excsamp = NULL)

## ----class_RAdata, eval=T-----------------------------------------------------
class(RAdata)

## ----RAdata, eval=T-----------------------------------------------------------
RAdata

## ----makeIC2, eval=T----------------------------------------------------------
samID = RAdata$extractVar("indID")$indID
head(samID)

## ----makeIC, eval=T-----------------------------------------------------------
progeny = samID[-c(1:4)]
mySpecies <- makeIC(RAobj = RAdata, samID = progeny,
                      filter = list(MAF = 0.05, MISS = 0.5,
                        BIN = 100, DEPTH = 5, PVALUE = 0.01, MAXDEPTH=500))

## ----2pt_rf-------------------------------------------------------------------
mySpecies$rf_2pt(nClust = 3)

## ----save, eval=FALSE---------------------------------------------------------
#  save(mySpecies, file="ICobject.Rdata")

## ----plotChr------------------------------------------------------------------
# plot 2-point rf matrix for BI and MI SNPs
mySpecies$plotChr(mat = "rf", lmai=0.5)

## ----plotChr2-----------------------------------------------------------------
# plot 2-point LOD matrix for BI, PI and MI SNPs
mySpecies$plotChr(mat = "LOD", lmai=0.5)

## ----createLG-----------------------------------------------------------------
mySpecies$createLG(LODthres = 10, nComp = 10)

## ----help_createLG, eval=F----------------------------------------------------
#  ?`$createLG`

## ----plotLG_mat---------------------------------------------------------------
mySpecies$plotLG(interactive = F, LG = NULL)

## ----mergeLG------------------------------------------------------------------
mySpecies$mergeLG(LG = c(3,4,5,7))
mySpecies$print(what = "LG")

## ----removeLG-----------------------------------------------------------------
mySpecies$removeLG(LG = 4)
mySpecies$print(what = "LG")

## ----removeSNP, eval=FALSE----------------------------------------------------
#  mySpecies$removeSNP(snps = c(36,27,32,28,131))

## ----createLG2----------------------------------------------------------------
mySpecies$createLG(LODthres = 8)

## ----orderLG, , fig.width=8, fig.height=6-------------------------------------
mySpecies$orderLG(mapfun = "haldane", weight = "LOD2", ndim = 30) 

## ----help_orderLG, eval=F-----------------------------------------------------
#  ?`$orderLG`

## ----computeMap---------------------------------------------------------------
mySpecies$computeMap(chrom = NULL, nThreads = 3, rfthres=0.1)

## ----computeMap2--------------------------------------------------------------
## remove SNP 173
mySpecies$removeSNP(c(124,133,89))
## Recompute linkage map
mySpecies$computeMap(chrom = c(1,2), nThreads = 3, rfthres=0.1)

## ----plotLM-------------------------------------------------------------------
mySpecies$plotLM(LG = NULL, fun = "haldane", col=c("cyan","lightblue","blue"))

## ----plotSyn, fig.width=7, fig.height=7---------------------------------------
mySpecies$plotSyn()

## ----writeLG, eval=F----------------------------------------------------------
#  mySpecies$writeLM(file = "mySpecies", inferGeno = TRUE)

## ----help_writeLM, eval=F-----------------------------------------------------
#  ?`$writeLM`

