% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addSNPs.R
\name{$addSNPs}
\alias{$addSNPs}
\title{BC and FS method: Add MI, PI and SI SNPs to existing linkage groups}
\usage{
BCobj$addSNPs(LODthres = 10, nComp = 30)
FSobj$addSNPs(LODthres = 10, nComp = 30)
}
\arguments{
\item{LODthres}{A positive numeric value specifying the LOD threshold used to add SNPs to the linkage groups.}

\item{nComp}{A positive integer value specifying how many SNPs in the linkage group to compute the average LOD score 
with the unmapped SNP.}
}
\description{
Method for adding maternal-informative (MI), paternal-informative (PI) and semi-informative (SI) SNPs 
to the paternal and maternal (pseudo-testcross) linkage groups.
}
\details{
Each MI, PI or SI SNP is mapped to the maternal and paternal linkage groups using the following
algorithm:
\enumerate{
\item For each pseudo-testcross linkage group, compute the average LOD score between the unmapped SNP and a 
specified number (\code{nComp}) of SNPs in the pseudo-testcross linkage group that have the highest LOD score 
with the unmapped SNP.
\item Map the unmapped SNP to the pseudo-testcross linkage group with the highest average LOD score if:
\itemize{
  \item The largest average LOD score is greater than the LOD threshold (\code{LODthres}).
  \item All the other average LOD scores associated with the other linkage groups are less than the LOD threshold.
  }
}
}
\examples{
## Simulate some sequencing data
set.seed(6745)
config <- list(list(sample(c(1,2,4), size=30, replace=TRUE)))
F1data <- simFS(0.01, config=config, meanDepth=10, nInd = 50)
## Compute 2-point recombination fractions
F1data$rf_2pt()
## create paternal and maternal linkage groups
F1data$createLG()

## Add more SNPs to the linkage groups
F1data$addSNPs()
}
\seealso{
\code{\link{BC}}, \code{\link{FS}}
}
\author{
Timothy P. Bilton
}
