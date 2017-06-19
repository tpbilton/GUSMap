\name{simFS}
\alias{simFS}
\title{Simulation of sequencing data of single full-sib family.}
\description{
Simulate genotypes and sequencing depth for sequencing data for a full-sib family population and write the data to various linkage mapping software formats.}
\usage{
  simFS <- function(rVec_f, rVec_m=rVec_f, config, nInd, nSnps, meanDepth, thres=NULL, NoDS=1,
                    formats=list(gusmap=T,onemap=F,lepmap=F,joinmap=F,crimap=F), rd_dist="Neg_Binom",
                    filename=NULL, seed1=1, seed2=1)
}
\arguments{
  \item{rVec_f, rVec_m}{vector of true paternal and maternal recombinaiton fractions. If only a single number is given,
  then recombination between adjacent loci are equal.}
  \item{config}{integer vector giveing the segregation type of each marker. 1=informative, 2=paternal segregating, 3=maternal segregating and 4=uniformative.}
  \item{nSnps}{number of loci}
  \item{nInd}{number of individuals}
  \item{meanDepth}{the mean depth for the read depth distribution.}
  \item{thres}{threshold value for which genotype calls with a depth less than the threshold are set to missing}
  \item{NoDS}{number of data sets to be simulated}
  \item{formats}{list where the element corresponding to the software name is TRUE will produce data in that format. Currently handles software for GusMap, CriMap, JoinMap, OneMap, LepMap}
  \item{rd_dist}{distribution for which the read depths are simulated from. Currently, only negative binomial ("Neg_Binom") and Poisson ("Pois") are implemented}
  \item{filename}{filename of the files generated from the simulation}
  \item{seed1}{random seed used for the simulation of the parental phase (or OPGP)}
  \item{seed2}{random seed used for the simulation of the data sets}
}
\details{
}
\value{}
\author{Timothy P. Bilton}
\examples{}
