\name{rf_est_FS}
\alias{rf_est_FS}
\title{Estimation of adjacent recombination fractions in full-sib families.}
\description{
Estimate the recombination fractions based on the hidden Markov model (HMM) for low coverage sequencing data in full-sib families.  
}
\usage{
rf_est_FS <- function(init_r=NULL, genon, depth, OPGP,
                      sexSpec=F, trace=F, noFam=1, \ldots)
}
\arguments{
  \item{init_r}{Vector of starting values for the recombination fractions}
  \item{genon}{List object with each element containing a genon matrix for each of the families.}
  \item{depth}{List object with each element containing a depth matrix for each of the families.}
  \item{OPGP}{Numeric vector of ordered parental genotype pairs (OPGPs) for each locus. See Bilton (2017) for a classification of each OPGP.}
  \item{sexSpec}{Logical value. If TRUE, sex spceific recombination fractions are estimated.}
  \item{trace}{Logical value. If TRUE, output from optim are printed.}
  \item{noFam}{Numeric value. Specifies the number of full-sib families used to estimate the recombination fractions.}
  \item{\ldots}{Additional arguments passed to optim}
}
\details{
The \code{genon} and \code{depth} matrices must have the rows representing the individuals and columns representing the SNPs. For entries of the \code{genon} matrix, 2=major homozygous, 1=heterozygous, 0=minor homozygous and NA missing genotype while the entries for the \code{depth} matrix must be a finite non-negative integer number. The number of families specified must match the number of genon and depths matrices given. For the vector of starting values, \code{init_r}, a single value can be given which sets all the starting values to be equal.

The likelihood calulations are scaled using forward recusion to advoid overflow issues and so can be implemented on a large numbers of loci. If the OPGP vector is unknown, then it can be inferred using \code{\link{infer_OPGP_FS}}.
}
\value{Function returns a list object. If non sex-specific recombination fractions are specified, the list contains,
\itemize{
\item{rf: }{Vector of recombination fraction estimates.}
\item{loglik: }{Value of the log-likelihood value at the recombination fraction estimates..}
}
Else, if sex-specific recombination fractions are desired, the list contains,
\itemize{
\item{rf_p: }{Vector of the paternal recombination fraction estimates.}
\item{rf_m: }{Vector of the maternal recombination fraction estimates.}
\item{loglik: }{Value of the log-likelihood value at the recombination fraction estimates.}
}
}
\references{
Bilton, T.P., Schofield, M.R., Black, M.A., Chagne D., Wilcox P., Dodds K.G. (2017). More accurate genetic maps using low coverage sequencing data in multiparental outcrossed populations. Unpublished manuscript.
}
\author{Timothy P. Bilton}
\seealso{
\code{\link{infer_OPGP_FS}}
}
\examples{
## simulate full sib family
config <- c(2,1,1,3,2,3,1,1,3,1,2,1)
F1data <- simFS(0.01, config=config, nInd=50, nSnps=12, meanDepth=5)

## Determine the parental phase
OPGP <- infer_OPGP_FS(F1data$genon, F1data$depth, config)

## Estimate the recombination fractions
rf_est_FS(genon = list(F1data$genon), depth = list(F1data$depth), OPGP = list(OPGP), noFam = 1)
}
  
