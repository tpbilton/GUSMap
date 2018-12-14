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

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix viterbi_fs_err(NumericVector rf, NumericVector ep, int nInd, int nSnps,
                             IntegerMatrix ref, IntegerMatrix alt, IntegerVector OPGP){
  
  // Compute the emission matrix
  NumericMatrix pAA(nInd, nSnps);
  NumericMatrix pAB(nInd, nSnps);
  NumericMatrix pBB(nInd, nSnps);
  int ind, snp, a, b;
  for(ind = 0; ind < nInd; ind++){
    for(snp = 0; snp < nSnps; snp++){
      a = ref(ind, snp);
      b = alt(ind, snp);
      pAA(ind, snp) = pow(1-ep[snp], a) * pow(ep[snp], b);
      pAB(ind, snp) = pow(0.5, a+b);
      pBB(ind, snp) = pow(1-ep[snp], b) * pow(ep[snp], a); 
    }
  }
  
  // run the viterbi algorithm
  int s1, s2, inferState;
  NumericMatrix eta(4, nSnps);
  NumericMatrix states(nInd, nSnps);
  NumericVector etaTemp(4);
  for(ind = 0; ind < nInd; ind++){
    // Compute forward probabilities at snp 1
    //sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      eta(s1, 0) = 0.25 * Qentry(OPGP[0], pAA(ind, 0), pAB(ind, 0), pBB(ind, 0), s1+1);
      //sum = sum + alphaDot[s1];
    }
    
    for(snp = 1; snp < nSnps; snp++){
      // compute the next forward probabilities for snp \ell
      for(s2 = 0; s2 < 4; s2++){
        for(s1 = 0; s1 < 4; s1++){
          etaTemp[s1] = eta(s1, snp) * Tmat(s1, s2, rf[snp-1]);
        }
        eta(s2, snp) = etaTemp[which_max(etaTemp)] * Qentry(OPGP[snp], pAA(ind, snp), pAB(ind, snp), pBB(ind, snp), s2+1);
      }
    }
    
/*    for(s1 = 0; s1 < 4; s1++){
      etaTemp[s1] = eta(s1, nSnps - 1);
    }
    inferState = which.max(etaTemp);
    states(ind, nSnps - 1) = inferState;
    // Now compute the most likely states in reverse
    for(snp = nSnps - 2; snp > -1; snp--){
      for(s1 = 0; s1 < 4; s1++){
        etaTemp[s1] = eta(s1, snp) * Tmat(s1, inferState, rf[snp-1]);
      }
      inferState = which.max(etaTemp);
      states(ind, snp) = inferState;
    }*/
  }
  return states;
}
  
  
  
