#include "Rcpp.h"

// Define the method signature

#ifdef __cplusplus
extern "C" {
#endif
  
#include "GUSMap.h"
  
#ifdef __cplusplus
}
#endif

Rcpp::NumericVector ll_fs_scaled_err_cpp(const Rcpp::NumericVector& r,
                                         const Rcpp::NumericVector& ep,
                                         const Rcpp::IntegerVector& ref,
                                         const Rcpp::IntegerVector& alt, 
                                         const Rcpp::NumericVector& bcoef_mat, 
                                         const Rcpp::NumericVector& Kab,
                                         const Rcpp::IntegerVector& OPGP,
                                         const Rcpp::NumericVector& nInd,
                                         const Rcpp::NumericVector& nSnps,
                                         const Rcpp::NumericVector& nThreads
                                         ) {
  
  SEXP ll = ll_fs_scaled_err_c(r, ep, ref, alt, bcoef_mat, Kab, OPGP, nInd, nSnps, nThreads);
  Rcpp::NumericVector result( ll );
  return result;
}