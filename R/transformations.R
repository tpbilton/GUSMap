### Transformation and inverse trnasformation for transforming the r.f's on the range [0,0.5] to [-Inf,Inf]
### Author: Timothy Bilton
### Date: 18/01/2017
### Edited: 23/03/2017

logit2 <- function(p) log(2*p/(1-2*p))

inv.logit2 <- function(logit_p) {
  p <- 1/(2*(1+1/exp(logit_p)))
  p.na <- is.na(p)
  if(sum(p.na)!=0)
    p[which(p.na)] <- 0.5
  return(p)
}
