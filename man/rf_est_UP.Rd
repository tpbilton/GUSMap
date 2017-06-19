\name{rf_est_FS_UP}
\alias{rf_est_FS_UP}
\title{To come}
\description{
  To come
}
\usage{
  rf_est_FS_UP <- function(genon, depth, config, trace=F, ...)
  }
\arguments{
  \item{genon}{matrix of genotypes. 2=major homozygous, 1=heterozygous, 0=minor homozygous}
  \item{depth}{matrix of read detphs}
  \item{config}{segregation type for each loci. 1=informative, 2=paternally segregating, 3=maternally segregating, 4=uniformative}
  \item{trace}{if TRUE, output from optim are printed}
  \item{...}{additional arguments passed  to optim}
}
\details{
}
\value{}
\author{Timothy P. Bilton}
\examples{}
  
  