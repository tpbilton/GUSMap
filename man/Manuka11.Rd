\name{Manuka11}
\alias{Manuka11}
\title{Manuka chromosome 11 SNPs}
\description{
Function for extracting the Manuka data used in the publication by Bilton et al. (2017).
}
\usage{
Manuka11 <- function()
}
\arguments{
}
\details{
The data consists fo 680 SNPs, genotyped using genotyping-by-sequencing methods. The original data is in RA format and the scripts from KGD (https://github.com/AgResearch/KGD) github repository are used to convert the data into a genon matrix and depth matrix as required by GusMap. 

}
\value{
Function outputs a list containing the following elements:
\itemize{
\item{genon: }{Numeric matrix of genotype calls. 2=major homozygous, 1=heterozygous, 0=minor homozygous and NA missing genotype}
\item{depth: }{Numeric matrix of sequencing depths.}
\item{indNames: }{Names of the individuals. The \emph{i}th entry corresponds to the \emph{i}th row of the \code{genon} and \code{depth} matrices.}
\item{chrom: }{Chromosome number of the SNPs. The \emph{i}th entry corresponds to the \emph{i}th column of the \code{genon} and \code{depth} matrices.}
\item{pos: }{Position number of the SNPs. The \emph{i}th entry corresponds to the \emph{i}th column of the \code{genon} and \code{depth} matrices.}
}
}
\author{Timothy P. Bilton}
\references{
Bilton, T.P., Schofield, M.R., Black, M.A., Chagne D., Wilcox P., Dodds K.G. (2017). Multilocus genetic linkage mapping in multiparental outcrossed populations using low coverage sequencing data. Unpublished manuscript.

Dodds K.G., McEwan J.C., Brauning R., Anderson R.M., van Stijn, T.C., Kristjansson, T., Clarke, S.M. (2015). Construction of relatedness matrices using genotyping-by-sequencing data. \emph{BMC Genomics}. 16:1047.
}
\examples{
## To convert the data into the required format for GusMap:
manuka <- Manuka11()

## SNPs removed due to mis-ordering
badSnps <- c(3:14,16,18:23,25:27,29,33:35,39:42,44,46:49,
  53,58:60,63,65:67,69:74,76,77,85,104,114:129,
  152,162,163,173:180,187:190,206:208,210,212,213,
  216,228,272,285,286,298,304,314,318,320,342,344,345,
  354,358,370,373,380,389,390,397,408,411,421,427,
  437,446,469,493,498:500,506,510,516,519,522,531,
  543,553,556,569,582,585,601,617,620,623,628,635,
  636:643,647:650,669:680)

## Work out the segregation types
apply(genon[1:4,],2,function(x){
if(all(x==2,na.rm=T)|all(x==0) & all(x==1))
})

## Compute the maps 
OPGP_rd6 <- infer_OPGP_FS(genon[,-badSnps], depth[,-badSnps], temp.config,
                          ndeps=rep(0.001,sum(temp.config==1)*2+sum(temp.config!=1)-2), reltol=1e-3)

## To access the RA directly
RAname <- system.file("extdata", "Manuka11.txt.gz", package="GusMap")
RAdata <- read.table(RAname, header=T)
}