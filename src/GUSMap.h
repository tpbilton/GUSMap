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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#ifndef _GUSMap
#define _GUSMap


SEXP ll_fs_scaled_err_c(SEXP r, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP OPGP, SEXP nInd, SEXP nSnps);
SEXP ll_fs_ss_scaled_err_c(SEXP r, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP OPGP, SEXP nInd, SEXP nSnps);
SEXP ll_fs_up_ss_scaled_err_c(SEXP r, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP config, SEXP nInd, SEXP nSnps);
SEXP EM_HMM(SEXP r, SEXP ep, SEXP depth_Ref, SEXP depth_Alt, SEXP OPGP, SEXP noFam, SEXP nInd, SEXP nSnps, SEXP sexSpec, SEXP seqError, SEXP para, SEXP ss_rf);
SEXP EM_HMM_UP(SEXP r, SEXP ep, SEXP depth_Ref, SEXP depth_Alt, SEXP config, SEXP noFam, SEXP nInd, SEXP nSnps, SEXP seqError, SEXP para, SEXP ss_rf);

#endif 
