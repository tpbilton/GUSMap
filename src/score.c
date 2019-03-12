/*
##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2019 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
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


// Derivative function for rf's
double der_rf(int s1, int s2, double rval){
  //double er;
  int sSum = s1 + s2*4;
  if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15)){
    return -2*rval*(1-rval)*(1-2*rval);
  }
  else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12)){
    return 2*rval*rval*(1-2*rval);
  }
  else{
    return rval*(1-2*rval)*(1-2*rval);
  }
}


// Derivative funcations for ep
double partial_der_ep(int geno, double ep, int a, int d){
  switch(geno){
  case 1:
      return binomial(a,d-a) *powl(ep,d-a) * powl(1-ep,a) * ((d-a)-d*ep);
  case 2:
    return 0;
  case 3:
      return binomial(a,d-a) * powl(ep,a) * powl(1-ep,d-a) * (a-d*ep);
  }
  return -1;
}


double der_ep(int OPGP, double ep, int a, int b, int elem){
  if((a == 0) & (b == 0))
    return 0;
  int d = a + b;
  switch(OPGP){
  case 1:
    if(elem == 1)
      return partial_der_ep(3, ep, a, d);
    else if ((elem == 2)|(elem == 3))  
      return 0;
    else if (elem == 4)
      return partial_der_ep(1, ep, a, d);
  case 2:
    if(elem == 3)
      return partial_der_ep(3, ep, a, d);
    else if ((elem == 1)|(elem == 4))
      return 0;
    else if (elem == 2)
      return partial_der_ep(1, ep, a, d);
  case 3:
    if(elem == 2) 
      return partial_der_ep(3, ep, a, d);
    else if ((elem == 1)|(elem == 4))
      return 0;
    else if (elem == 3)
      return partial_der_ep(1, ep, a, d);
  case 4:
    if(elem == 4) 
      return partial_der_ep(3, ep, a, d);
    else if ((elem == 2)|(elem == 3))
      return 0;
    else if (elem == 1)
      return partial_der_ep(1, ep, a, d);
  case 5:
    if ((elem == 1)|(elem == 2))
      return 0;
    else if ((elem == 3)|(elem == 4))
      return partial_der_ep(1, ep, a, d);
  case 6:
    if ((elem == 1)|(elem == 2))
      return partial_der_ep(1, ep, a, d);
    else if ((elem == 3)|(elem == 4))
      return 0;
  case 7:
    if ((elem == 1)|(elem == 2))
      return partial_der_ep(3, ep, a, d);
    else if ((elem == 3)|(elem == 4))
      return 0;
  case 8:
    if ((elem == 1)|(elem == 2))
      return 0;
    else if ((elem == 3)|(elem == 4))
      return partial_der_ep(3, ep, a, d);
  case 9:
    if ((elem == 1)|(elem == 3))
      return 0;
    else if ((elem == 2)|(elem == 4))
      return partial_der_ep(1, ep, a, d);
  case 10:
    if ((elem == 1)|(elem == 3))
      return partial_der_ep(1, ep, a, d);
    else if ((elem == 2)|(elem == 4))
      return 0;
  case 11:
    if ((elem == 1)|(elem == 3))
      return partial_der_ep(3, ep, a, d);
    else if ((elem == 2)|(elem == 4))
      return 0;
  case 12:
    if ((elem == 1)|(elem == 3))
      return 0;
    else if ((elem == 2)|(elem == 4))
      return partial_der_ep(3, ep, a, d);
  case 13:
    return partial_der_ep(1, ep, a, d);
  case 14:
    return 0;
  case 15:
    return 0;
  case 16:
    return partial_der_ep(3, ep, a, d);
  } // end of Switch
  return -999;
}

// rf's are equal and single error parameter
SEXP score_fs_scaled_multi_err_c(SEXP r, SEXP ep, SEXP ref, SEXP alt, SEXP bcoef_mat, SEXP Kab, 
                                 SEXP OPGP, SEXP nInd, SEXP nSnps, SEXP nThreads){
  // Initialize variables
  int ind, snp, snp_der, nInd_c, nSnps_c, nThreads_c, *pOPGP, *pref, *palt, maxThreads;
  double *pscore, *pr, *pKab, *pep, *pbcoef_mat;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = INTEGER(OPGP);
  pref = INTEGER(ref);
  palt = INTEGER(alt);
  pKab = REAL(Kab);
  pbcoef_mat = REAL(bcoef_mat);
  pr = REAL(r);
  pep = REAL(ep);
  // Define the output variable
  SEXP score;
  PROTECT(score = allocVector(REALSXP, 2*nSnps_c-1));
  pscore = REAL(score);
  //SEXP pout = PROTECT(allocVector(VECSXP, 3));
  
  // set up number of threads
  nThreads_c = asInteger(nThreads);
  maxThreads = omp_get_max_threads();
  if (nThreads_c <= 0) {
    // if nThreads is set to zero then use everything
    nThreads_c = maxThreads;
  }
  else if (nThreads_c > maxThreads) {
    // don't allow more threads than the maximum available
    nThreads_c = maxThreads;
  }
  
  // define the density values for the emission probs
  int index;
  double pKaa[nSnps_c * nInd_c];
  double pKbb[nSnps_c * nInd_c];
//  #pragma omp parallel for num_threads(nThreads_c) private(index, snp)
  for(ind = 0; ind < nInd_c; ind++){
    for(snp = 0; snp < nSnps_c; snp++){
      index = ind + snp*nInd_c;
      pKaa[index] = pbcoef_mat[index] * pow(1.0 - pep[snp], pref[index]) * pow(pep[snp], palt[index]);
      pKbb[index] = pbcoef_mat[index] * pow(1.0 - pep[snp], palt[index]) * pow(pep[snp], pref[index]);
    }
  }
  
  double llval = 0, score_c[2*nSnps_c-1];
  for(snp = 0; snp < 2*nSnps_c-1; snp++){
    score_c[snp] = 0;
  }
  
  // Now compute the likelihood and score function
  #pragma omp parallel for reduction(+:llval) num_threads(nThreads_c) private(snp, snp_der)
  for(ind = 0; ind < nInd_c; ind++){
    int s1, s2;
    double phi[4][2*nSnps_c-1], phi_prev[4][2*nSnps_c-1];
    double alphaTilde[4], alphaDot[4], sum, sum_der, w_new, w_prev;
    double delta;
    
    // Compute forward probabilities at snp 1
    sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      alphaDot[s1] = 0.25 * Qentry(pOPGP[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1);
      sum = sum + alphaDot[s1];
      alphaTilde[s1] = alphaDot[s1];
      // Compute the derivative for ep
      phi_prev[s1][nSnps_c-1] = 0.25 * der_ep(pOPGP[0], pep[0], pref[ind], palt[ind], s1+1);
    }
    
    // add contribution to likelihood
    w_prev = sum;
    
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps_c; snp++){
      // compute the next forward probabilities for snp \ell
      w_new = 0;
      for(s2 = 0; s2 < 4; s2++){
        sum = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum = sum + Tmat(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        delta = Qentry(pOPGP[snp], pKaa[ind + nInd_c*snp], pKab[ind + nInd_c*snp], pKbb[ind + nInd_c*snp], s2+1);
        alphaDot[s2] = sum * delta/w_prev;
        // add contribution to new weight
        w_new = w_new + alphaDot[s2];
        // Compute the derivatives
        //rf's
        sum_der = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum_der = sum_der + der_rf(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        phi[s2][snp-1] = sum_der * delta * 1/w_prev;
        for(snp_der = 0; snp_der < snp-1; snp_der++){
          sum_der = 0;
          for(s1 = 0; s1 < 4; s1++){
            sum_der = sum_der + phi_prev[s1][snp_der] * Tmat(s1, s2, pr[snp-1]);
          }
          phi[s2][snp_der] = sum_der * delta * 1/w_prev;
        }
        //sequencing error parameter
        sum_der = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum_der = sum_der + alphaTilde[s1] * der_ep(pOPGP[snp], pep[snp], pref[ind + nInd_c*snp], palt[ind + nInd_c*snp], s2+1) * Tmat(s1, s2, pr[snp-1]);
        }
        phi[s2][nSnps_c-1 + snp] = sum_der * 1/w_prev;
        for(snp_der = 0; snp_der < snp; snp_der++){
          sum_der = 0;
          for(s1 = 0; s1 < 4; s1++){  
            sum_der = sum_der + phi_prev[s1][snp_der + nSnps_c-1] * Tmat(s1, s2, pr[snp-1]);
          }
          phi[s2][snp_der + nSnps_c-1] = sum_der * delta * 1/w_prev;
        }
      }
      // Add contribution to the likelihood
      llval = llval + log(w_prev);
      
      // Scale the forward probability vector
      for(s2 = 0; s2 < 4; s2++){
        alphaTilde[s2] = alphaDot[s2];
        // update the derivative vectors
        for(snp_der = 0; snp_der < 2*nSnps_c-1; snp_der++){
          phi_prev[s2][snp_der] = phi[s2][snp_der];
        }
      }
      w_prev = w_new;
    }
    llval = llval + log(w_prev);
    // add contributions to the score vector
    for(snp_der = 0; snp_der < 2*nSnps_c-1; snp_der++){
      sum_der = 0;
      for(s2 = 0; s2 < 4; s2++)
        sum_der = sum_der + phi[s2][snp_der]/w_prev;
      #pragma omp atomic
      score_c[snp_der] += sum_der;
    }
  }
  // Compute the score for each parameter
  for(snp_der=0; snp_der < 2*nSnps_c-1; snp_der++){
    pscore[snp_der] = score_c[snp_der];
  }
  
  // Clean up and return likelihood value
  UNPROTECT(1);
  return score;
}
    
    
    
// rf's are equal and single error parameter
SEXP score_fs_scaled_err_c(SEXP r, SEXP ep, SEXP ref, SEXP alt, SEXP bcoef_mat, SEXP Kab, 
        SEXP OPGP, SEXP nInd, SEXP nSnps, SEXP nThreads){
  // Initialize variables
  int ind, snp, snp_der, nInd_c, nSnps_c, nThreads_c, *pOPGP, *pref, *palt, maxThreads;
  double *pscore, *pr, *pKab, ep_c, *pbcoef_mat;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = INTEGER(OPGP);
  pref = INTEGER(ref);
  palt = INTEGER(alt);
  pKab = REAL(Kab);
  pbcoef_mat = REAL(bcoef_mat);
  pr = REAL(r);
  ep_c = REAL(ep)[0];
  // Define the output variable
  SEXP score;
  PROTECT(score = allocVector(REALSXP, nSnps_c));
  pscore = REAL(score);
  //SEXP pout = PROTECT(allocVector(VECSXP, 3));

  // set up number of threads
  nThreads_c = asInteger(nThreads);
  maxThreads = omp_get_max_threads();
  if (nThreads_c <= 0) {
    // if nThreads is set to zero then use everything
    nThreads_c = maxThreads;
  }
  else if (nThreads_c > maxThreads) {
    // don't allow more threads than the maximum available
    nThreads_c = maxThreads;
  }
  
  // define the density values for the emission probs
  long size = (long) nSnps_c * nInd_c;
  double pKaa[size];
  double pKbb[size];
  #pragma omp parallel for num_threads(nThreads_c)
  for (long i = 0; i < size; i++) {
    pKaa[i] = pbcoef_mat[i] * pow(1.0 - ep_c, pref[i]) * pow(ep_c, palt[i]);
    pKbb[i] = pbcoef_mat[i] * pow(1.0 - ep_c, palt[i]) * pow(ep_c, pref[i]);
  }
  
  double llval = 0, score_c[nSnps_c];
  for(snp = 0; snp < nSnps_c; snp++){
    score_c[snp] = 0;
  }

  // Now compute the likelihood and score function
  #pragma omp parallel for reduction(+:llval) num_threads(nThreads_c) \
                           private(snp, snp_der)
  for(ind = 0; ind < nInd_c; ind++){
    int s1, s2;
    double phi[4][nSnps_c], phi_prev[4][nSnps_c];
    double alphaTilde[4], alphaDot[4], sum, sum_der, w_new, w_prev;
    double delta;

    // Compute forward probabilities at snp 1
    sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      alphaDot[s1] = 0.25 * Qentry(pOPGP[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1);
      sum = sum + alphaDot[s1];
      alphaTilde[s1] = alphaDot[s1];
      // Compute the derivative for ep
      phi_prev[s1][nSnps_c-1] = 0.25 * der_ep(pOPGP[0], ep_c, pref[ind], palt[ind], s1+1);
    }

    // add contribution to likelihood
    w_prev = sum;

    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps_c; snp++){
      // compute the next forward probabilities for snp \ell
      w_new = 0;
      for(s2 = 0; s2 < 4; s2++){
        sum = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum = sum + Tmat(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        delta = Qentry(pOPGP[snp], pKaa[ind + nInd_c*snp], pKab[ind + nInd_c*snp], pKbb[ind + nInd_c*snp], s2+1);
        alphaDot[s2] = sum * delta/w_prev;
        // add contribution to new weight
        w_new = w_new + alphaDot[s2];
        // Compute the derivatives
        //rf's
        sum_der = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum_der = sum_der + der_rf(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        phi[s2][snp-1] = sum_der * delta * 1/w_prev;
        for(snp_der = 0; snp_der < snp-1; snp_der++){
          sum_der = 0;
          for(s1 = 0; s1 < 4; s1++){
            sum_der = sum_der + phi_prev[s1][snp_der] * Tmat(s1, s2, pr[snp-1]);
          }
          phi[s2][snp_der] = sum_der * delta * 1/w_prev;
        }
        //sequencing error parameter
        sum_der = 0;
        for(s1 = 0; s1 < 4; s1++){  
          sum_der = sum_der + ((phi_prev[s1][nSnps_c-1] * delta + alphaTilde[s1] * 
            der_ep(pOPGP[snp], ep_c, pref[ind + nInd_c*snp], palt[ind + nInd_c*snp], s2+1)) * Tmat(s1, s2, pr[snp-1]));
        }
        phi[s2][nSnps_c-1] = sum_der/w_prev;
      }
      // Add contribution to the likelihood
      llval = llval + log(w_prev);

      // Scale the forward probability vector
      for(s2 = 0; s2 < 4; s2++){
        alphaTilde[s2] = alphaDot[s2];
        // update the derivative vectors
        for(snp_der = 0; snp_der < nSnps_c; snp_der++){
          phi_prev[s2][snp_der] = phi[s2][snp_der];
        }
      }
      w_prev = w_new;
    }
    llval = llval + log(w_prev);
    //Rprintf("Likelihood value: %f\n", llval);
    // add contributions to the score vector
    for(snp_der = 0; snp_der < snp-1; snp_der++){
      sum_der = 0;
      for(s2 = 0; s2 < 4; s2++)
        sum_der = sum_der + phi[s2][snp_der]/w_prev;
      #pragma omp atomic
      score_c[snp_der] += sum_der;
    }
    for(snp_der = nSnps_c-1; snp_der < nSnps_c; snp_der++){
      sum_der = 0;
      for(s2 = 0; s2 < 4; s2++)
        sum_der = sum_der + phi[s2][snp_der]/w_prev;
      #pragma omp atomic
      score_c[snp_der] += sum_der;
    }
  }
  // Compute the score for each parameter
  for(snp_der=0; snp_der < nSnps_c; snp_der++){
    pscore[snp_der] = score_c[snp_der];
  }
  
  // Clean up and return likelihood value
  //Rprintf("Likelihood value: %f\n", llval);
  UNPROTECT(1);
  return score;
}

// rf's are equal and no error parameter
SEXP score_fs_scaled_c(SEXP r, SEXP ep, SEXP ref, SEXP alt, SEXP bcoef_mat, SEXP Kab,
                       SEXP OPGP, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int s1, s2, ind, snp, snp_der, nInd_c, nSnps_c, *pOPGP, *palt, *pref;
  double *pscore, *pr, *pbcoef_mat, *pKab, ep_c, delta;
  double alphaTilde[4], alphaDot[4], sum, sum_der, w_new, w_prev;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pref = INTEGER(ref);
  palt = INTEGER(alt);
  pOPGP = INTEGER(OPGP);
  pKab = REAL(Kab);
  pbcoef_mat = REAL(bcoef_mat);
  pr = REAL(r);
  ep_c = REAL(ep)[0];
  // Define the output variable
  SEXP score;
  PROTECT(score = allocVector(REALSXP, nSnps_c-1));
  pscore = REAL(score);
  //SEXP pout = PROTECT(allocVector(VECSXP, 3));
  double llval = 0, phi[4][nSnps_c-1], phi_prev[4][nSnps_c-1], score_c[nSnps_c-1];
  for(snp = 0; snp < nSnps_c-1; snp++){
    score_c[snp] = 0;
  }
  
  // define the density values for the emission probs
  long size = (long) nSnps_c * nInd_c;
  double pKaa[size];
  double pKbb[size];
  // #pragma omp parallel for
  for (long i = 0; i < size; i++) {
    pKaa[i] = pbcoef_mat[i] * pow(1.0 - ep_c, pref[i]) * pow(ep_c, palt[i]);
    pKbb[i] = pbcoef_mat[i] * pow(1.0 - ep_c, palt[i]) * pow(ep_c, pref[i]);
  }
  
  // Now compute the likelihood and score function
  for(ind = 0; ind < nInd_c; ind++){
    // Compute forward probabilities at snp 1
    sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      alphaDot[s1] = 0.25 * Qentry(pOPGP[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1);
      sum = sum + alphaDot[s1];
      alphaTilde[s1] = alphaDot[s1];
    }

    // add contribution to likelihood
    w_prev = sum;
    
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps_c; snp++){
      // compute the next forward probabilities for snp \ell
      w_new = 0;
      for(s2 = 0; s2 < 4; s2++){
        sum = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum = sum + Tmat(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        delta = Qentry(pOPGP[snp], pKaa[ind + nInd_c*snp], pKab[ind + nInd_c*snp], pKbb[ind + nInd_c*snp], s2+1);
        alphaDot[s2] = sum * delta/w_prev;
        // add contribution to new weight
        w_new = w_new + alphaDot[s2];
        // Compute the derivatives
        //rf's
        sum_der = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum_der = sum_der + der_rf(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        phi[s2][snp-1] = sum_der * delta * 1/w_prev;
        for(snp_der = 0; snp_der < snp-1; snp_der++){
          sum_der = 0;
          for(s1 = 0; s1 < 4; s1++){
            sum_der = sum_der + phi_prev[s1][snp_der] * Tmat(s1, s2, pr[snp-1]);
          }
          phi[s2][snp_der] = sum_der * delta * 1/w_prev;
        }
      }
      // Add contribution to the likelihood
      llval = llval + log(w_prev);
      
      // Scale the forward probability vector
      for(s2 = 0; s2 < 4; s2++){
        alphaTilde[s2] = alphaDot[s2];
        // update the derivative vectors
        for(snp_der = 0; snp_der < nSnps_c-1; snp_der++){
          phi_prev[s2][snp_der] = phi[s2][snp_der];
        }
      }
      w_prev = w_new;
    }
    llval = llval + log(w_prev);
    //Rprintf("Likelihood value: %f\n", llval);
    // add contributions to the score vector
    for(snp_der = 0; snp_der < snp-1; snp_der++){
      sum_der = 0;
      for(s2 = 0; s2 < 4; s2++)
        sum_der = sum_der + phi[s2][snp_der]/w_prev;
      sum_der = sum_der + score_c[snp_der];
      score_c[snp_der] = sum_der;
    }
  }
  // Compute the score for each parameter
  for(snp_der=0; snp_der < nSnps_c-1; snp_der++){
    pscore[snp_der] = score_c[snp_der];
  }
  
  // Clean up and return likelihood value
  //Rprintf("Likelihood value: %f\n", llval);
  UNPROTECT(1);
  return score;
}

