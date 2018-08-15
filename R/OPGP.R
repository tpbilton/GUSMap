##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping
# Copyright 2017-2018 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################
#### Functions for converting phase information to OPGPs
#### Author: Timothy P. Bilton
#### Date: 06/02/16


## Function for deteriming the parental phase of a full-sib family
#' Inference of the OPGPs (or parental phase) for a single full-sub family.
#' 
#' Infers the OPGPs for all loci of a single full-sib family.
#' 
#' The \code{ref} and \code{alt} matrices must have the rows
#' representing the individuals and columns representing the SNPs. The entries
#' of these two matrices must be a finite non-negative integer number. The
#' value for the \code{config} vector are 1=both-informative (ABxAB),
#' 2=paternal-informative A (ABxAA), 3=paternal-informative B (ABxBB),
#' 4=maternal-informative (AAxAB), 5=maternal-informative (BBxAB), 6 =
#' uninformative (AAxAA), 7 = uninformative (AAxBB), 8 = uninformative (BBxAA)
#' and 9 = uninformative (BBxBB).
#' 
#' The methodology used to infer the OPGP is found in the section "Inferring
#' OPGPs" of the manuscript by Bilton et al. (2017). Note that inference is
#' made on a single full-sib family at a time.
#' 
#' Two different optimzation procedures are available, which are the EM algorithm and optim.
#' To control the parameters to these procedures, addition arguments can be passed to the function.
#' The arguments which have an effect are dependent on the optimization procedure.
#' \itemize{
#' \item EM: Only two arguments currently have an effect. 'reltol' specifies 
#' the maximum difference between the likelihood value of successive iterations
#' before the algorithm terminates. 'maxit' specifies the maximum number of iterations
#' used in the algorithm
#' \item optim: The extra arguments are passed directly to optim. Those see what 
#' arguments are valid, visit the help page fro optim using '?optim'.
#' }
#' 
#' Note: There can be issues at times in finding the maximum likelihood
#' estimate (MLE) of the likelihood for the implementation which uses optim.
#' Changing the \code{ndeps} argument that is passed to \code{optim()} sometimes helps.
#' Alternatively, one can just use the EM approach (which is the default).
#' 
#' @param ref Numeric matrix of allele counts for the reference allele.
#' @param alt Numeric matrix of allele counts for the alternate allele.
#' @param config Numeric vector of the segregation types for all the loci.
#' @param ep Numeric value of the starting value for the sequencing error
#' parameter. Default is NULL which means that the sequencing parameter is
#' fixed at zero.
#' @param method A character string specifying whether to use the EM algorithm
#' or optim to perform the optimzation.
#' @param \ldots Additional arguments passed to optimization procedure. See details for more information.
#' @return Function returns a vector of the inferred OPGP values. These values correspond to those 
#' given in Table 1 of Bilton et al. (2017).
#' @author Timothy P. Bilton
#' @references Bilton, T.P., Schofield, M.R., Black, M.A., Chagné, D., Wilcox,
#' P.L., Dodds K.G. (2017). Accounting for errors in low coverage high-throughput
#' sequencing data when constructing genetic maps using biparental outcrossed
#' populations. Unpublished Manuscript.
#' @examples
#' 
#' ## simulate full-sib family
#' config <- c(2,1,1,4,2,4,1,1,4,1,2,1)
#' F1data <- simFS(0.01, config=config, nInd=50, meanDepth=5)
#' 
#' ## Determine the parental phase
#' OPGP <- infer_OPGP_FS(F1data$ref, F1data$alt, config)
#' # Same and the true OPGP value 
#' F1data$OPGP
#' 
#' ## if there are issues with solving the likelihood (see comment in details)
#' ## Note: Only if using optim method
#' OPGP <- infer_OPGP_FS(F1data$ref, F1data$alt, config, method = "optim",
#'                       ndeps=rep(1e-1,sum(config==1)*2+sum(config!=1)-1))
#' 
#' @export infer_OPGP_FS

infer_OPGP_FS <- function(ref, alt, config, ep=0.001, method="optim", nThreads=0, ...){
  
  if(!is.matrix(ref) || !is.matrix(alt))
    stop("The read counts inputs are not matrix objects")
  if( (!is.null(ep) & !is.numeric(ep)) )
    stop("Starting values for the error parameter needs to be a single numeric value in the interval (0,1) or a NULL object")
  if( !is.numeric(config) || !is.vector(config) || any(!(config == round(config))) || any(config < 1) || any(config > 9) )
    stop("Segregation information needs to be an integer vector equal to the number of SNPs with entires from 1 to 9")
  if(length(method) != 1 || !(method %in% c("EM","optim")) )
    stop("Optimzation method specified is invalid. Please select one of 'EM' or 'optim'")
  
  nSnps <- ncol(ref); nInd <- nrow(ref)
  
  if(nInd*nSnps > 25000)          # if data set is too large, there are memory issues with R for EM algorithm
    method = "optim"
  
  ## Index the informative loci
  Isnps <- which(!(config %in% 6:9))
  
  ## solve the likelihood for when phase is not known
  MLEs <- rf_est_FS_UP(ref[,Isnps], alt[,Isnps], config[Isnps], ep=ep, method=method, nThreads=nThreads, ...)
  
  parHap <- matrix("A",nrow=4,ncol=nSnps)
  parHap[cbind(c(rep(3,sum(config %in% c(3,7,9))),rep(4,sum(config %in% c(3,7,9)))),which(config %in% c(3,7,9)))] <- "B"
  parHap[cbind(c(rep(1,sum(config %in% c(5,8,9))),rep(2,sum(config %in% c(5,8,9)))),which(config %in% c(5,8,9)))] <- "B"
  ## Determine whether the phase between adjacent paternal (or maternal) informative SNPs is in coupling
  coupling_pat <- MLEs[[1]] < 0.5
  coupling_mat <- MLEs[[2]] < 0.5
  
  # Work out the indices of the r.f. parameters of each sex
  ps_snps <- which(config %in% c(1,2,3))
  ms_snps <- which(config %in% c(1,4,5))
  npar <- c(length(ps_snps)-1,length(ms_snps)-1)
  
  s_pat <- numeric(npar[[1]]+1)
  s_pat[1] <- 2
  
  s_mat <- numeric(npar[[2]]+1)
  s_mat[1] <- 2
  
  for(i in 1:npar[[1]]){
    s_pat[i+1] <- ifelse(coupling_pat[i],s_pat[i],s_pat[i]%%2+1)
  }
  for(i in 1:npar[[2]]){
    s_mat[i+1] <- ifelse(coupling_mat[i],s_mat[i],s_mat[i]%%2+1)
  }
  
  parHap[cbind(s_pat,ps_snps)] <- "B"
  parHap[cbind(s_mat+2,ms_snps)] <- "B"
  
  return(parHapToOPGP(parHap))
}

## Function for converting parental haplotypes to OPGPs
parHapToOPGP <- function(parHap, major="A", minor="B"){
  return (apply(parHap,2,function(x){
    if (all(x==major)){
      return(13) 
    } else if (x[1]==x[2] & x[3]==x[4] & x[1]==major & x[3]==minor){
      return(14)
    } else if (x[1]==x[2] & x[3]==x[4] & x[3]==major & x[1]==minor){
      return(15)
    } else if (all(x==minor)){
      return(16)
    } else if(sum(x==major)==2){
      return((x[1]==minor)+(x[3]==minor)*2+1)
    } else if (all(x[3:4]==major)){
      return((x[1]==minor)+5)
    } else if (all(x[1:2]==major)){
      return((x[3]==minor)+9)
    } else if (all(x[3:4]==minor)){
      return((x[1]==minor)+7)
    } else if (all(x[1:2]==minor)){
      return((x[3]==minor)+11)
    } 
  }) 
)}



