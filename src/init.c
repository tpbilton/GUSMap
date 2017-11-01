#include "GUSMap.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static const R_CallMethodDef callMethods[] = {
  {"ll_fs_scaled_err_c",       (DL_FUNC) &ll_fs_scaled_err_c,		7},
  {"ll_fs_ss_scaled_err_c",    (DL_FUNC) &ll_fs_ss_scaled_err_c,	7},
  {"ll_fs_up_ss_scaled_err_c", (DL_FUNC) &ll_fs_up_ss_scaled_err_c,	7},
  {NULL,		   NULL,				0}
};

void R_init_GUSMap(DllInfo *info){
  
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);

  R_useDynamicSymbols(info, TRUE);

  R_RegisterCCallable("GUSMap","ll_fs_scaled_err_c", 		(DL_FUNC) &ll_fs_scaled_err_c);
  R_RegisterCCallable("GUSMap","ll_fs_ss_scaled_err_c", 	(DL_FUNC) &ll_fs_ss_scaled_err_c);
  R_RegisterCCallable("GUSMap","ll_fs_up_ss_scaled_err_c",      (DL_FUNC) &ll_fs_up_ss_scaled_err_c);
}
