\name{rf_est_FS}
\alias{rf_est_FS}
\title{Estimation of adjacent recombination fractions in full-sib families.}
\description{
Estimate the recombination fractions based on the hidden Markov model (HMM) for low coverage sequencing data in full-sib families.  
}
\usage{
rf_est_FS(init_r=NULL, epsilon=NULL, depth_Ref, depth_Alt, OPGP,
                      sexSpec=F, trace=F, noFam=1, \ldots)
}
\arguments{
  \item{init_r}{Vector of starting values for the recombination fractions}
  \item{epsilon}{Numeric value of the starting value for the sequencing error parameter. Default is NULL which means that the sequencing parameter is fixed at zero.}
  \item{depth_Ref}{List object with each element containing a matrix of allele counts for the reference allele for each family.}
  \item{depth_Alt}{List object with each element containing a matrix of allele counts for the alternate allele for each family.}
  \item{OPGP}{List object with each element containing a numeric vector of ordered parental genotype pairs (OPGPs) for each family. See Bilton (2017) for a classification of each OPGP.}
  \item{sexSpec}{Logical value. If TRUE, sex specific recombination fractions are estimated.}
  \item{trace}{Logical value. If TRUE, output from \code{optim()} are printed.}
  \item{noFam}{Numeric value. Specifies the number of full-sib families used to estimate the recombination fractions.}
  \item{\ldots}{Additional arguments passed to \code{optim()}.}
}
\details{
The \code{depth_Ref} and \code{depth_Alt} matrices must have the rows representing the individuals and columns representing the SNPs. The entries of these two matrices must be a finite non-negative integer number. The number of families specified must match the number of genon (genotype calls) and depth matrices given. For the vector of starting values, \code{init_r}, a single value can be given which sets all the starting values to be equal.

The likelihood calculations are scaled using forward recursion to avoid overflow issues and so can be implemented on a large numbers of loci. If the OPGP vector is unknown for the families, then it can be inferred using \code{\link{infer_OPGP_FS}}.
}
\value{Function returns a list object. If non sex-specific recombination fractions are specified, the list contains,
\itemize{
\item rf: Vector of recombination fraction estimates.
\item epsilon: Estimate of the sequencing error parameter.
\item loglik: The log-likelihood value at the maximum likelihood estimates.
}
Else, if sex-specific recombination fractions are desired, the list contains,
\itemize{
\item rf_p: Vector of the paternal recombination fraction estimates.
\item rf_m: Vector of the maternal recombination fraction estimates.
\item epsilon: Estimate of the sequencing error parameter.
\item loglik: The log-likelihood value at the maximum likelihood estimates.
}
}
\references{
Bilton, T.P., Schofield, M.R., Black, M.A., Chagne D., Wilcox P., Dodds K.G. (2017). Accounting for errors in low coverage sequencing data when constructing genetic maps using biparental outcrossed populations. Unpublished Manuscript.
}
\author{Timothy P. Bilton}
\seealso{
\code{\link{infer_OPGP_FS}}
}
\examples{
### Case 1: Single family
## simulate full sib family
config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
F1data <- simFS(0.01, config=config, nInd=50, meanDepth=5)

## Determine the parental phase
OPGP <- infer_OPGP_FS(F1data$depth_Ref, F1data$depth_Alt, config)

## Estimate the recombination fractions
rf_est_FS(depth_Ref = list(F1data$depth_Ref), depth_Alt = list(F1data$depth_Alt), OPGP = list(OPGP), noFam = 1)

########################
### Case 2: Two families
config_1 <- c(1,1,3,1,4,4,1,6,2,4,6,1)
config_2 <- c(2,4,1,6,2,1,1,2,5,6,1,1)

Fam1 <- simFS(0.01, config=config_1, nInd=50, meanDepth=5)
Fam2 <- simFS(0.01, config=config_2, nInd=50, meanDepth=5)

## Determine the parental phase in each family
OPGP_1 <- infer_OPGP_FS(Fam1$depth_Ref, Fam1$depth_Alt, config_1)
OPGP_2 <- infer_OPGP_FS(Fam2$depth_Ref, Fam2$depth_Alt, config_2)

## Estimate the recombination fractions
rf_est_FS(depth_Ref = list(Fam1$depth_Ref,Fam2$depth_Ref),
          depth_Alt = list(Fam1$depth_Alt,Fam2$depth_Alt), OPGP = list(OPGP_1,OPGP_2), noFam = 2)
}
  
