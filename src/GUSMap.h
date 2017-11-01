#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifndef _GUSMap
#define _GUSMap


SEXP ll_fs_scaled_err_c(SEXP r, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP OPGP, SEXP nInd, SEXP nSnps);
SEXP ll_fs_ss_scaled_err_c(SEXP r, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP OPGP, SEXP nInd, SEXP nSnps);
SEXP ll_fs_up_ss_scaled_err_c(SEXP r, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP config, SEXP nInd, SEXP nSnps);

#endif 
