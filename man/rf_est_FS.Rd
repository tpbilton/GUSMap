\name{rf_est_FS}
\alias{rf_est_FS}
\title{Estimation of adjacent recombination fractions using a HMM}
\description{
To come
}
\usage{
  rf_est_FS <- function(init_r=NULL, genon, depth, OPGP, sexSpec=F, trace=F, scaled=F, hessian=F, noFam=1, ...){
}
\arguments{
  \item{init_r}{starting values for the recombination fractions.}
  \item{genon}{matrix of genotypes. 2=major homozygous, 1=heterozygous, 0=minor homozygous}
  \item{depth}{matrix of read detphs}
  \item{OPGP}{Ordered parental genotype pair (OPGP) for each locus}
  \item{sexSpec}{if TRUE, sex spceific recombination fractions are estimated}
  \item{trace}{if TRUE, output from optim are printed.}
  \item{noFam}{number of families that are going to be used in the estimation procedure}
  \item{...}{additional arguments passed  to optim}
}
\details{
}
\value{}
\author{Timothy P. Bilton}
\examples{}
  
