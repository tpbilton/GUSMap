\name{rf_est_FS_UP}
\alias{rf_est_FS_UP}
\title{Estimation of sex-specific adjacent recombination fractions in full-sib families when the parental phase is unknown.}
\description{
Estimates the recombination fractions for the hidden Markov model (HMM) for low coverage sequencing data in full-sib families when the paretnal phase is unknown. 
}
\usage{
rf_est_FS_UP <- function(genon, depth, config, trace=F, \ldots)
}
\arguments{
  \item{genon}{Numeric matrix of genotype calls.}
  \item{depth}{Numeric matrix of sequencing depths.}
  \item{config}{Numeric vector of the segregation types for all the loci.}
  \item{trace}{Logical value. If TRUE, output from optim are printed.}
  \item{\ldots}{Additional arguments passed to optim.}
}
\details{
The \code{genon} and \code{depth} matrices must have the rows representing the individuals and columns representing the SNPs. For entries of the \code{genon} matrix, 2=major homozygous, 1=heterozygous, 0=minor homozygous and NA missing genotype while the entries for \code{depth} matrix must be a finite non-negative integer number.  The value for the \code{config} vector are 1=informative, 2=paternally segregating, 3=maternally segregating and 4=uniformative.

The main difference of this function from \link{\code{rf_est_FS}} is that the recombination fractions are constrainted to the interval [0,1] and are foced to be sex-specific. Primiarly, this function is included to infer phase and is used by the function \link{\code{infer_OPGP_FS}} although it can be used to estimate 2 point recombination fraction estimates between partially segregating loci.
}
\value{
The function returns a list with elements
estimated recombation fraction estimates and }
\author{Timothy P. Bilton}
\seealso{
\code{\link{infer_OPGP_FS}}
}
\examples{
## simulate full sib family
config <- c(2,2)
F1data <- simFS(0.01, config=config, nInd=50, nSnps=2, meanDepth=5)

## Estimate 2-point r.f.'s
rf_est_FS_UP(F1data$genon, F1data$depth, config)
}
  
  