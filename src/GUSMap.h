#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifndef _GusMap
#define _GusMap


SEXP ll_fs_c(SEXP r, SEXP genon, SEXP depth, SEXP OPGP, SEXP nInd, SEXP nSnps);
SEXP ll_fs_ss_c(SEXP r, SEXP genon, SEXP depth, SEXP OPGP, SEXP nInd, SEXP nSnps);
SEXP ll_fs_up_ss_c(SEXP r, SEXP genon, SEXP depth, SEXP config, SEXP nInd, SEXP nSnps);
SEXP ll_fs_scaled_c(SEXP r, SEXP genon, SEXP depth, SEXP OPGP, SEXP nInd, SEXP nSnps);
SEXP ll_fs_ss_scaled_c(SEXP r, SEXP genon, SEXP depth, SEXP OPGP, SEXP nInd, SEXP nSnps);
SEXP ll_fs_up_ss_scaled_c(SEXP r, SEXP genon, SEXP depth, SEXP config, SEXP nInd, SEXP nSnps);

#endif 
