/*
##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
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
 */

#include "GUSMap.h"
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static const R_CallMethodDef callMethods[] = {
  {"ll_fs_scaled_err_c",       (DL_FUNC) &ll_fs_scaled_err_c,		    8},
  {"ll_fs_ss_scaled_err_c",    (DL_FUNC) &ll_fs_ss_scaled_err_c,	  7},
  {"ll_fs_up_ss_scaled_err_c", (DL_FUNC) &ll_fs_up_ss_scaled_err_c,	7},
  {"score_fs_scaled_c",        (DL_FUNC) &score_fs_scaled_c,        7},
  {"score_fs_scaled_err_c",    (DL_FUNC) &score_fs_scaled_err_c,    11},
  {"EM_HMM",                   (DL_FUNC) &EM_HMM,                   13},
  {"EM_HMM_UP",                (DL_FUNC) &EM_HMM_UP,                12},
  {NULL,		                   NULL,	                              0}
};

void R_init_GUSMap(DllInfo *info){
  
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);

  R_useDynamicSymbols(info, TRUE);

  R_RegisterCCallable("GUSMap","ll_fs_scaled_err_c", 		        (DL_FUNC) &ll_fs_scaled_err_c);
  R_RegisterCCallable("GUSMap","ll_fs_ss_scaled_err_c", 	      (DL_FUNC) &ll_fs_ss_scaled_err_c);
  R_RegisterCCallable("GUSMap","ll_fs_up_ss_scaled_err_c",      (DL_FUNC) &ll_fs_up_ss_scaled_err_c);
  R_RegisterCCallable("GUSMap","score_fs_scaled_c",             (DL_FUNC) &score_fs_scaled_c);
  R_RegisterCCallable("GUSMap","score_fs_scaled_err_c",         (DL_FUNC) &score_fs_scaled_err_c);
  R_RegisterCCallable("GUSMap","EM_HMM",                        (DL_FUNC) &EM_HMM);
  R_RegisterCCallable("GUSMap","EM_HMM_UP",                     (DL_FUNC) &EM_HMM_UP);
}
