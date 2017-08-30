\name{infer_OPGP_FS}
\alias{infer_OPGP_FS}
\title{Inference of the OPGPs (or parental phase) for full-sub families.}
\description{
Infers the OPGPs for all loci of a single full-sib family.
}
\usage{
infer_OPGP_FS <- function(genon, depth, config, \ldots)
}
\arguments{
  \item{genon}{Numeric matrix of genotype calls.}
  \item{depth}{Numeric matrix of sequencing depths.}
  \item{config}{Numeric vector of the segregation types for all the loci. 1=doubly informative, 2=paternal segregating, 3=maternal segregating and 4=uniformative.}
  \item{\ldots}{Additional arguments passed to optim.}
}
\details{
The \code{genon} and \code{depth} matrices must have the rows representing the individuals and columns representing the SNPs. For entries of the \code{genon} matrix, 2=major homozygous, 1=heterozygous, 0=minor homozygous and NA missing genotype while the entries for the \code{depth} matrix must be a finite non-negative integer number. The value for the \code{config} vector are 1=informative, 2=paternally segregating, 3=maternally segregating and 4=uniformative.

Note: There can be issues at times in finding the maximum likelihood estimate (MLE) of the likelihood. Changing the \code{ndeps} argument that is passed to optim tends to help.
}
\value{Function returns a vector of the inferred OPGP values}
\author{Timothy P. Bilton}
\seealso{
\code{\link{rf_est_UP}}
}
\examples{
## simulate full-sib family
config <- c(2,1,1,3,2,3,1,1,3,1,2,1)
F1data <- simFS(0.01, config=config, nInd=50, nSnps=12, meanDepth=5)

## Determine the parental phase
OPGP <- infer_OPGP_FS(F1data$genon, F1data$depth, config)

## if there are issues with solving the likelihood (see comment in details)
OPGP <- infer_OPGP_FS(F1data$genon, F1data$depth, config,
                      ndeps=rep(1e-1,sum(config==1)*2+sum(config!=1)-2))
}
  
  