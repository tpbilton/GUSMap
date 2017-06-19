#include "GusMap.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static const R_CallMethodDef callMethods[] = {
  {"ll_fs_c",		   (DL_FUNC) &ll_fs_c,			6},
  {"ll_fs_ss_c",	   (DL_FUNC) &ll_fs_ss_c,		6},
  {"ll_fs_up_ss_c",	   (DL_FUNC) &ll_fs_up_ss_c,		6},
  {"ll_fs_scaled_c", 	   (DL_FUNC) &ll_fs_scaled_c,		6},
  {"ll_fs_ss_scaled_c",    (DL_FUNC) &ll_fs_ss_scaled_c,	6},
  {"ll_fs_up_ss_scaled_c", (DL_FUNC) &ll_fs_up_ss_scaled_c,	6},
  {NULL,		   NULL,				0}
};

void R_init_GBSanalysis(DllInfo *info){
  
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);

  R_useDynamicSymbols(info, TRUE);

  R_RegisterCCallable("GusMap","ll_fs_c", 			(DL_FUNC) &ll_fs_c);
  R_RegisterCCallable("GusMap","ll_fs_ss_c", 		(DL_FUNC) &ll_fs_ss_c);
  R_RegisterCCallable("GusMap","ll_fs_up_ss_c", 		(DL_FUNC) &ll_fs_up_ss_c);
  R_RegisterCCallable("GusMap","ll_fs_scaled_c", 		(DL_FUNC) &ll_fs_scaled_c);
  R_RegisterCCallable("GusMap","ll_fs_ss_scaled_c", 	(DL_FUNC) &ll_fs_ss_scaled_c);
  R_RegisterCCallable("GusMap","ll_fs_up_ss_scaled_c",  (DL_FUNC) &ll_fs_up_ss_scaled_c);
}
