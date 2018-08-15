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
#include <Rmath.h>
#include <math.h>
#include "probFun.h"
#ifdef _OPENMP
    #include <omp.h>
#else
    inline int omp_get_max_threads() { return 1; }
#endif

//////////// likelihood functions for multipoint likelihood in full sib-families using GBS data /////////////////////
// Input variables for likelihoods
//  - r: recombination fraction values
//  - genon: matrix of genotype calls. Format in C is a 1-D vector
//  - depth: matrix of depth values. Format in C is a 1-D vector
//  - OPGP: The OPGP of all the SNPs.
//  - nInd: Number of individuals.
//  - nSnps: Number of SNPs.
//  - config: Parental genotype configurations (or segregation type)
//            1 = Informative, 2 = paternal segregating, 3 = maternal segregating


//// Scaled versions of the likelihood to deal with overflow issues


//// likelihood 3:
// Not sex-specific (assumed equal) 
// r.f constrainted to range [0,1/2].
// OPGP's (or phase) are assumed to be known
// Include error parameters
SEXP ll_fs_scaled_err_c(SEXP r, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP OPGP, SEXP nInd, SEXP nSnps, SEXP nThreads){
  // Initialize variables
  int s1, s2, ind, snp, nInd_c, nSnps_c, nThreads_c, *pOPGP;
  double *pll, *pr, *pKaa, *pKab, *pKbb;
  double alphaTilde[4], alphaDot[4], sum, w_new;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  nThreads_c = asInteger(nThreads);
  // Define the pointers to the other input R variables
  pOPGP = INTEGER(OPGP);
  pKaa = REAL(Kaa);
  pKab = REAL(Kab);
  pKbb = REAL(Kbb);
  pr = REAL(r);  
  // Define the output variable
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  pll = REAL(ll);
  double llval = 0;

  // if nThreads is set to zero then use everything
  if (nThreads_c <= 0) {
    nThreads_c = omp_get_max_threads();
  }
  
  // Now compute the likelihood
  #pragma omp parallel for reduction(+:llval) num_threads(nThreads_c) \
                           private(sum, s1, alphaDot, alphaTilde, snp, s2, w_new)
  for(ind = 0; ind < nInd_c; ind++){
    // Compute forward probabilities at snp 1
    sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      //Rprintf("Q value :%.6f at snp %i in ind %i\n", Qentry(pOPGP[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1, delta_c), 0, ind);
      alphaDot[s1] = 0.25 * Qentry(pOPGP[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1);
      sum = sum + alphaDot[s1];
    }
    // Scale forward probabilities
    for(s1 = 0; s1 < 4; s1++){
      alphaTilde[s1] = alphaDot[s1]/sum;
    }
    //Rprintf("New weight :%.6f at snp %i\n", sum, 1);
    
    // add contribution to likelihood
    //w_logcumsum = log(sum);
    llval = llval + log(sum);
    
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps_c; snp++){
      // compute the next forward probabilities for snp \ell
      for(s2 = 0; s2 < 4; s2++){
        sum = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum = sum + Tmat(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        //Rprintf("Q value :%.6f at snp %i in ind %i\n", Qentry(pOPGP[snp], pKaa[ind + nInd_c*snp], pKab[ind + nInd_c*snp], pKbb[ind + nInd_c*snp], s2+1, delta_c), snp, ind);
        alphaDot[s2] = Qentry(pOPGP[snp], pKaa[ind + nInd_c*snp], pKab[ind + nInd_c*snp], pKbb[ind + nInd_c*snp], s2+1) * sum;
      }
      // Compute the weight for snp \ell
      w_new = 0;
      for(s2 = 0; s2 < 4; s2++){
        w_new = w_new + alphaDot[s2];
      }
      // Add contribution to the likelihood
      llval = llval + log(w_new);
      //w_logcumsum = w_logcumsum + log(w_new);
      // Scale the forward probability vector
      for(s2 = 0; s2 < 4; s2++){
        alphaTilde[s2] = alphaDot[s2]/w_new;
      }
    }
  }

  pll[0] = -1*llval;
  // Clean up and return likelihood value
  UNPROTECT(1);
  return ll;
}

//// likelihood 2:
// sex-specific
// r.f constrainted to range [0,1/2].
// OPGP's (or phase) are assumed to be known
SEXP ll_fs_ss_scaled_err_c(SEXP r, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP OPGP, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int s1, s2, ind, snp, nInd_c, nSnps_c, *pOPGP;
  double *pll, *pr, *pKaa, *pKab, *pKbb;
  double alphaTilde[4], alphaDot[4], sum, w_new;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = INTEGER(OPGP);
  pKaa = REAL(Kaa);
  pKab = REAL(Kab);
  pKbb = REAL(Kbb);
  pr = REAL(r);  
  // Define the output variable
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  pll = REAL(ll);
  double llval = 0;
  
  // Now compute the likelihood
  for(ind = 0; ind < nInd_c; ind++){
    // Compute forward probabilities at snp 1
    sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      alphaDot[s1] = 0.25 * Qentry(pOPGP[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1);
      sum = sum + alphaDot[s1];
    }
    // Scale forward probabilities
    for(s1 = 0; s1 < 4; s1++){
      alphaTilde[s1] = alphaDot[s1]/sum;
    }
    //Rprintf("New weight :%.6f at snp %i\n", sum, 1);
    
    // add contribution to likelihood
    //w_logcumsum = log(sum);
    llval = llval + log(sum);
    
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps_c; snp++){
      // compute the next forward probabilities for snp \ell
      for(s2 = 0; s2 < 4; s2++){
        sum = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum = sum + Tmat_ss(s1, s2, pr[snp-1], pr[snp-1+nSnps_c-1]) * alphaTilde[s1];
        }
        alphaDot[s2] = Qentry(pOPGP[snp], pKaa[ind + nInd_c*snp], pKab[ind + nInd_c*snp], pKbb[ind + nInd_c*snp], s2+1) * sum;
      }
      // Compute the weight for snp \ell
      w_new = 0;
      for(s2 = 0; s2 < 4; s2++){
        w_new = w_new + alphaDot[s2];
      }
      // Add contribution to the likelihood
      llval = llval + log(w_new);
      //w_logcumsum = w_logcumsum + log(w_new);
      // Scale the forward probability vector
      for(s2 = 0; s2 < 4; s2++){
        alphaTilde[s2] = alphaDot[s2]/w_new;
      }
    }
  }
  
  pll[0] = -1*llval;
  // Clean up and return likelihood value
  UNPROTECT(1);
  return ll;
}


//// likelihood 3: For sex-specific unphased data
// sex-specific
// r.f constrainted to range [0,1].
// OPGP's (or phase) are not known
SEXP ll_fs_up_ss_scaled_err_c(SEXP r, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP config, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int s1, s2, ind, snp, nInd_c, nSnps_c, *pconfig;
  double *pll, *pr, *pKaa, *pKab, *pKbb;
  double alphaTilde[4], alphaDot[4], sum, w_new;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pconfig = INTEGER(config);
  pKaa = REAL(Kaa);
  pKab = REAL(Kab);
  pKbb = REAL(Kbb);
  pr = REAL(r);  
  // Define the output variable
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  pll = REAL(ll);
  double llval = 0;
  
  // Now compute the likelihood
  for(ind = 0; ind < nInd_c; ind++){
    // Compute forward probabilities at snp 1
    sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      //Rprintf("Q value :%.6f at snp %i in ind %i\n", Qentry_up(pconfig[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1, delta_c), 0, ind);
      alphaDot[s1] = 0.25 * Qentry_up(pconfig[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1);
      sum = sum + alphaDot[s1];
    }
    // Scale forward probabilities
    for(s1 = 0; s1 < 4; s1++){
      alphaTilde[s1] = alphaDot[s1]/sum;
    }
    //Rprintf("New weight :%.6f at snp %i\n", sum, 1);
    
    // add contribution to likelihood
    //w_logcumsum = log(sum);
    llval = llval + log(sum);
    
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps_c; snp++){
      // compute the next forward probabilities for snp \ell
      for(s2 = 0; s2 < 4; s2++){
        sum = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum = sum + Tmat_ss(s1, s2, pr[snp-1], pr[snp-1+nSnps_c-1]) * alphaTilde[s1];
        }
        //Rprintf("Q value :%.6f at snp %i in ind %i\n", Qentry_up(pconfig[0], pKaa[ind + nInd_c*snp], pKab[ind + nInd_c*snp], pKbb[ind + nInd_c*snp], s2+1, delta_c), snp, ind);
        alphaDot[s2] = Qentry_up(pconfig[snp], pKaa[ind + nInd_c*snp], pKab[ind + nInd_c*snp], pKbb[ind + nInd_c*snp], s2+1) * sum;
      }
      // Compute the weight for snp \ell
      w_new = 0;
      for(s2 = 0; s2 < 4; s2++){
        w_new = w_new + alphaDot[s2];
      }
      // Add contribution to the likelihood
      llval = llval + log(w_new);
      //w_logcumsum = w_logcumsum + log(w_new);
      // Scale the forward probability vector
      for(s2 = 0; s2 < 4; s2++){
        alphaTilde[s2] = alphaDot[s2]/w_new;
      }
    }
  }
  
  pll[0] = -1*llval;
  // Clean up and return likelihood value
  UNPROTECT(1);
  return ll;
}
