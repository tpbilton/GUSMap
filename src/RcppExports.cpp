// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// viterbi_fs_err
Rcpp::NumericMatrix viterbi_fs_err(Rcpp::NumericVector rf, Rcpp::NumericVector ep, int nInd, int nSnps, Rcpp::IntegerMatrix ref, Rcpp::IntegerMatrix alt, Rcpp::IntegerVector OPGP);
RcppExport SEXP _GUSMap_viterbi_fs_err(SEXP rfSEXP, SEXP epSEXP, SEXP nIndSEXP, SEXP nSnpsSEXP, SEXP refSEXP, SEXP altSEXP, SEXP OPGPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rf(rfSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type ep(epSEXP);
    Rcpp::traits::input_parameter< int >::type nInd(nIndSEXP);
    Rcpp::traits::input_parameter< int >::type nSnps(nSnpsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type ref(refSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type alt(altSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type OPGP(OPGPSEXP);
    rcpp_result_gen = Rcpp::wrap(viterbi_fs_err(rf, ep, nInd, nSnps, ref, alt, OPGP));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP EM_HMM(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP EM_HMM_multierr(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP EM_HMM_UP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP ll_fs_scaled_err_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP ll_fs_ss_scaled_err_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP ll_fs_up_ss_scaled_err_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP score_fs_scaled_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP score_fs_scaled_err_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP score_fs_scaled_multi_err_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP score_fs_ss_scaled_err_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP score_fs_ss_scaled_multi_err_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_GUSMap_viterbi_fs_err", (DL_FUNC) &_GUSMap_viterbi_fs_err, 7},
    {"EM_HMM",                         (DL_FUNC) &EM_HMM,                         13},
    {"EM_HMM_multierr",                (DL_FUNC) &EM_HMM_multierr,                13},
    {"EM_HMM_UP",                      (DL_FUNC) &EM_HMM_UP,                      12},
    {"ll_fs_scaled_err_c",             (DL_FUNC) &ll_fs_scaled_err_c,             10},
    {"ll_fs_ss_scaled_err_c",          (DL_FUNC) &ll_fs_ss_scaled_err_c,           7},
    {"ll_fs_up_ss_scaled_err_c",       (DL_FUNC) &ll_fs_up_ss_scaled_err_c,       10},
    {"score_fs_scaled_c",              (DL_FUNC) &score_fs_scaled_c,               9},
    {"score_fs_scaled_err_c",          (DL_FUNC) &score_fs_scaled_err_c,          10},
    {"score_fs_scaled_multi_err_c",    (DL_FUNC) &score_fs_scaled_multi_err_c,    10},
    {"score_fs_ss_scaled_err_c",       (DL_FUNC) &score_fs_ss_scaled_err_c,       11},
    {"score_fs_ss_scaled_multi_err_c", (DL_FUNC) &score_fs_ss_scaled_multi_err_c, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_GUSMap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
