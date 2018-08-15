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
#include <stdio.h>
#include <limits.h>
#include "probFun.h"

#ifdef _OPENMP
    #include <omp.h>
#else
    inline int omp_get_max_threads() { return 1; }
#endif

int Tcount(int s1, int s2){
  int sSum = s1 + s2*4;
  if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15))
    return 0;
  else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12))
    return 2;
  else
    return 1;
}


/*
double long binomial(int a, int b){
  double long res;
  res = expl(lgammal(a+b+1) - lgammal(b+1) - lgammal(a+1));
  return res;
}
*/

// Function for computing the emission probabilities given the true genotypes (parallel)
double computeProb(int nInd, int nSnps, double pAA[nInd][nSnps], double pBB[nInd][nSnps], double bin_coef[nInd][nSnps],
                        double epsilon, int *pref, int *palt, int nThreads_c){
  int ind;
  #pragma omp parallel for num_threads(nThreads_c)
  for(ind = 0; ind < nInd; ind++){
    int snp;
    for(snp = 0; snp < nSnps; snp++){
      int indx = ind + nInd * snp;
      if( (pref[indx] + palt[indx]) == 0){
        pAA[ind][snp] = 1;
        pBB[ind][snp] = 1;
      }
      else{
        pAA[ind][snp] = bin_coef[ind][snp] * powl(1 - epsilon, pref[indx]) * powl(epsilon, palt[indx]);
        pBB[ind][snp] = bin_coef[ind][snp] * powl(epsilon, pref[indx]) * powl(1 - epsilon, palt[indx]);
      }
    }
  }
  return -1;
}

// Function for extracting entries of the emission probability matrix
// when the OPGPs are known
int Iindx(int OPGP, int elem){
  switch(OPGP){
  case 1:
    if(elem == 1)
      return 2;
    else if ((elem == 2)|(elem == 3))  
      return 0;
    else if (elem == 4)
      return 1;
  case 2:
    if(elem == 3)
      return 2;
    else if ((elem == 1)|(elem == 4))
      return 0;
    else if (elem == 2)
      return 1;
  case 3:
    if(elem == 2) 
      return 2;
    else if ((elem == 1)|(elem == 4))
      return 0;
    else if (elem == 3)
      return 1;
  case 4:
    if(elem == 4) 
      return 2;
    else if ((elem == 2)|(elem == 3))
      return 0;
    else if (elem == 1)
      return 1;
  case 5:
    if ((elem == 1)|(elem == 2))
      return 0;
    else if ((elem == 3)|(elem == 4))
      return 1;
  case 6:
    if ((elem == 1)|(elem == 2))
      return 1;
    else if ((elem == 3)|(elem == 4))
      return 0;
  case 7:
    if ((elem == 1)|(elem == 2))
      return 2;
    else if ((elem == 3)|(elem == 4))
      return 0;
  case 8:
    if ((elem == 1)|(elem == 2))
      return 0;
    else if ((elem == 3)|(elem == 4))
      return 2;
  case 9:
    if ((elem == 1)|(elem == 3))
      return 0;
    else if ((elem == 2)|(elem == 4))
      return 1;
  case 10:
    if ((elem == 1)|(elem == 3))
      return 1;
    else if ((elem == 2)|(elem == 4))
      return 0;
  case 11:
    if ((elem == 1)|(elem == 3))
      return 2;
    else if ((elem == 2)|(elem == 4))
      return 0;
  case 12:
    if ((elem == 1)|(elem == 3))
      return 0;
    else if ((elem == 2)|(elem == 4))
      return 2;
  case 13:
    return 1;
  case 14:
    return 0;
  case 15:
    return 0;
  case 16:
    return 2;
  } // end of Switch
  return -1;
}




SEXP EM_HMM(SEXP r, SEXP ep, SEXP ref, SEXP alt, SEXP OPGP, SEXP noFam, SEXP nInd, SEXP nSnps,
            SEXP sexSpec, SEXP seqError, SEXP para, SEXP ss_rf, SEXP nThreads){
  // Initialize variables
  int s1, s2, fam, ind, snp, g, iter, nIter, parent, noFam_c, nSnps_c, sexSpec_c, seqError_c, nThreads_c;
  double sum, sumA, sumB, a, b, delta;

  // set up number of threads
  nThreads_c = asInteger(nThreads);
  if (nThreads_c <= 0) {
    // if nThreads is set to zero then use everything
    nThreads_c = omp_get_max_threads();
  }

  // Copy values of R input into C objects
  nSnps_c = INTEGER(nSnps)[0];
  double r_c[(nSnps_c-1)*2], ep_c;
  for(snp = 0; snp < nSnps_c - 1; snp++){
    for(parent = 0; parent < 2; parent++){
      r_c[snp+parent*(nSnps_c-1)] = REAL(r)[snp+parent*(nSnps_c-1)];
    }
  }
  ep_c = REAL(ep)[0];
  sexSpec_c = INTEGER(sexSpec)[0];
  seqError_c = INTEGER(seqError)[0];
  noFam_c = INTEGER(noFam)[0];
  int nInd_c[noFam_c], nTotal, indSum[noFam_c];
  nTotal = 0;
  for(fam = 0; fam < noFam_c; fam++){
    nInd_c[fam] = INTEGER(nInd)[fam];
    indSum[fam] = nTotal;
    nTotal = nTotal + nInd_c[fam];
    //indSum[fam] = nTotal - nInd_c[fam];
  }
  nIter = REAL(para)[0];
  delta = REAL(para)[1];
  // Initialize some more variables
  double alphaTilde[4][nTotal][nSnps_c], alphaDot[4];
  double betaTilde[4][nTotal][nSnps_c], betaDot[4];
  double w_new, log_w[nTotal][nSnps_c];
  double uProb[4][nTotal][nSnps_c], vProb[4][4][nTotal][nSnps_c-1];
  // Define the pointers to the other input R variables
  int *pOPGP, *pref, *palt;
  pOPGP = INTEGER(OPGP);
  pref = INTEGER(ref);
  palt = INTEGER(alt);
  int *pss_rf;
  pss_rf = INTEGER(ss_rf);
  // if sex-specific rf
  if(sexSpec_c){
    for(snp = 0; snp < nSnps_c-1; snp++){
      if(pss_rf[snp]==0){
        r_c[snp] = 0;
      }
      if(pss_rf[snp+nSnps_c-1]==0){
        r_c[snp + nSnps_c-1] = 0;
      }
    }
  } 
  // Probability matrics
  double pAA[nTotal][nSnps_c];
  double pAB[nTotal][nSnps_c];
  double pBB[nTotal][nSnps_c];
  double bin_coef[nTotal][nSnps_c];
  int Iaa[4][nTotal][nSnps_c], Ibb[4][nTotal][nSnps_c];
  //double (*pbin_coef)[nSnps_c] = &bin_coef;
  // Compute the probs for the heterozygous calls
  for(fam = 0; fam < noFam_c; fam++){
    for(ind = 0; ind < nInd_c[fam]; ind++){
      int indx = ind + indSum[fam];
      for(snp = 0; snp < nSnps_c; snp++){
        bin_coef[indx][snp] = 1; //binomial(pref[indx + nTotal*snp], palt[indx + nTotal*snp]);
        pAB[indx][snp] = bin_coef[indx][snp] * powl(0.5,pref[indx + nTotal*snp] + palt[indx + nTotal*snp]);
        for(s1 = 0; s1 < 4; s1++){
          g = Iindx(pOPGP[snp*noFam_c + fam], s1 + 1);
          if( g == 1 ){
            Iaa[s1][indx][snp] = 1;
            Ibb[s1][indx][snp] = 0;
          }
          else if( g == 2 ){
            Ibb[s1][indx][snp] = 1;
            Iaa[s1][indx][snp] = 0;
          }
          else{
            Iaa[s1][indx][snp] = 0;
            Ibb[s1][indx][snp] = 0;
          }
        }
      }
    }
  }
  // Define the output variable
  double *prout, *pepout, *pllout;
  SEXP rout = PROTECT(allocVector(REALSXP, 2*(nSnps_c-1)));
  SEXP epout = PROTECT(allocVector(REALSXP, 1));
  SEXP llout = PROTECT(allocVector(REALSXP, 1));
  prout = REAL(rout);
  pepout = REAL(epout);
  pllout = REAL(llout);
  SEXP pout = PROTECT(allocVector(VECSXP, 3));
  double llval = 0, prellval = 0;

  /////// Start algorithm
  iter = 0;
  while( (iter < 2) || ((iter < nIter) & ((llval - prellval) > delta))){
    iter = iter + 1;
    prellval = llval;
    llval = 0;
    
    // unpdate the probabilites for pAA and pBB given the parameter values and data.
    computeProb(nTotal, nSnps_c, pAA, pBB, bin_coef, ep_c, pref, palt, nThreads_c);
  
    // Compute the forward and backward probabilities for each individual
    for(fam = 0; fam < noFam_c; fam++){
      #pragma omp parallel for reduction(+:llval) num_threads(nThreads_c) \
                               private(sum, s1, s2, alphaDot, snp, w_new, betaDot)
      for(ind = 0; ind < nInd_c[fam]; ind++){
        int indx = ind + indSum[fam];
        //Rprintf("indSum %i\n", indSum[fam]);
        //Rprintf("indx %i\n", indx);
        //Rprintf("bin_coef %f at snp %i: ind %i:\n",bin_coef[indx][0], 0, indx);
        //Rprintf("pAA %f at snp %i: ind %i:\n",pAA[indx][0], 0, indx);
        //Rprintf("pAB %f at snp %i: ind %i:\n",pAB[indx][0], 0, indx);
        //Rprintf("pBB %f at snp %i: ind %i:\n",pBB[indx][0], 0, indx);
        /////////////////////////////////////////
        // Compute forward probabilities at snp 1
        sum = 0;
        //Rprintf("OPGP :%.6f at snp %i in ind %i\n", pOPGP[fam], 0, ind);
        for(s1 = 0; s1 < 4; s1++){
          //Rprintf("pAA %.16f, pAB %.16f pBB %.16f at snp %i in ind %i\n", pAA[indx][0], pAB[indx][0], pBB[indx][0], 0, ind);
          //Rprintf("Q value :%.22f at snp %i in ind %i\n", Qentry(pOPGP[fam], pAA[indx][0], pAB[indx][0], pBB[indx][0], s1+1), 0, ind);
          alphaDot[s1] = 0.25 * Qentry(pOPGP[fam], pAA[indx][0], pAB[indx][0], pBB[indx][0], s1+1);
          sum = sum + alphaDot[s1];
        }
        //Rprintf("Q value :%.22f at snp %i in ind %i\n", sum, 0, ind);
        // Scale forward probabilities
        for(s1 = 0; s1 < 4; s1++){
          alphaTilde[s1][indx][0] = alphaDot[s1]/sum;
        }
        //Rprintf("New weight :%.6f at snp %i\n", sum, 1);
      
        // add contribution to likelihood
        log_w[indx][0] = log(sum);
        llval = llval + log(sum);
        //w_logcumsum = log(sum);
        //llval = llval + w_logcumsum;
        // llval = llval - log(sum);
        //Rprintf("llvalue :%.6f at snp %i\n", llval, 0);
        
        // iterate over the remaining SNPs
        for(snp = 1; snp < nSnps_c; snp++){
          //Rprintf("snp %i and ind %i\n", snp, ind);
          //Rprintf("snp %i\n", snp );
          //Rprintf("nSnps %i\n", nSnps_c );
          // compute the next forward probabilities for snp j
          //Rprintf("pAA %f at snp %i: ind %i:\n",pAA[indx][snp], snp, ind);
          //Rprintf("pAB %f at snp %i: ind %i:\n",pAB[indx][snp], snp, ind);
          //Rprintf("pBB %f at snp %i: ind %i:\n",pBB[indx][snp], snp, ind);
          for(s2 = 0; s2 < 4; s2++){
            sum = 0;
            for(s1 = 0; s1 < 4; s1++){
              sum = sum + Tmat_ss(s1, s2, r_c[snp-1], r_c[snp-1+nSnps_c-1]) * alphaTilde[s1][indx][snp-1];
              
              //Rprintf("Tmat %.16f for r_c %.12f and %.12f and state %i and %i\n", Tmat_ss(s1, s2, r_c[snp-1], r_c[snp-1+nSnps_c-1]), r_c[snp-1], r_c[snp-1+nSnps_c-1], s1, s2);
            }
            //Rprintf("sum :%.16f at state %i, snp %i and ind %i\n", sum, s2, snp, ind);
            //Rprintf("Q value :%.22f at snp %i in ind %i\n", Qentry(pOPGP[snp*noFam_c + fam], pAA[indx][snp], pAB[indx][snp], pBB[indx][snp], s2+1), snp, ind);
            //Rprintf("pAA %.16f, pAB %.16f pBB %.16f at snp %i in ind %i\n", pAA[indx][snp], pAB[indx][snp], pBB[indx][snp], snp, indx);
            //Rprintf("Q value :%.22f at snp %i in ind %i\n", Qentry(pOPGP[snp*noFam_c + fam], pAA[indx][snp], pAB[indx][snp], pBB[indx][snp], s2+1), snp, indx);
            alphaDot[s2] = Qentry(pOPGP[snp*noFam_c + fam], pAA[indx][snp], pAB[indx][snp], pBB[indx][snp], s2+1) * sum;
            //Rprintf("alphaDot :%.22f at state %i, snp %i and ind %i\n", alphaDot[s2], s2, snp, ind);
          }
          
          // Compute the weight for snp \ell
          w_new = 0;
          for(s2 = 0; s2 < 4; s2++){
            w_new = w_new + alphaDot[s2];
          }
          // Add contribution to the likelihood
          //llval = llval + log(w_new) + w_logcumsum;
          //w_logcumsum = w_logcumsum + log(w_new);
          log_w[indx][snp] = log(w_new);
          // Scale the forward probability vector
          for(s2 = 0; s2 < 4; s2++){
            alphaTilde[s2][indx][snp] = alphaDot[s2]/w_new;
            //Rprintf("alphaTilde :%.22f at state %i, snp %i and ind %i\n", alphaTilde[s2][indx][snp], s2, snp, ind);
          }
          // Add contribution to the likelihood
          llval = llval + log(w_new);
          //Rprintf("llvalue :%.8f at iter %i\n", llval, indx);
        }
        //////////////////////
        // Compute the backward probabilities
        for(s1 = 0; s1 < 4; s1++){
          betaTilde[s1][indx][nSnps_c-1] = 1/exp(log_w[indx][nSnps_c-1]);
          //sum = log(betaTilde[s1][indx][nSnps_c-1]) - (log(1) -  log_w[indx][nSnps_c-1]);
          //Rprintf("log(beta) :%.22f at state %i, snp %i and ind %i\n", exp(sum), s1, nSnps_c-1, ind);
        }
        // iterate over the remaining SNPs
        sum = 0;
        for(snp = nSnps_c-2; snp > -1; snp--){
          for(s1 = 0; s1 < 4; s1++){
            sum = 0;
            for(s2 = 0; s2 < 4; s2++){
              sum = sum + Qentry(pOPGP[(snp+1)*noFam_c + fam], pAA[indx][snp+1], pAB[indx][snp+1], pBB[indx][snp+1], s2+1) * 
                Tmat_ss(s1, s2, r_c[snp], r_c[snp+nSnps_c-1]) * betaTilde[s2][indx][snp+1];
            }
            betaDot[s1] = sum;
          }
          // Scale the backward probability vector
          for(s1 = 0; s1 < 4; s1++){
            betaTilde[s1][indx][snp] = betaDot[s1]/exp(log_w[indx][snp]);
            //Rprintf("betaTilde :%.6f at state %i, snp %i and ind %i\n", betaTilde[s1][indx][snp], s1, snp, ind);
          }
        }
        ///////// E-step:
        for(snp = 0; snp < nSnps_c - 1; snp++){
          for(s1 = 0; s1 < 4; s1++){
            uProb[s1][indx][snp] = (alphaTilde[s1][indx][snp] * betaTilde[s1][indx][snp])/exp(-log_w[indx][snp]);
            for(s2 = 0; s2 < 4; s2++){
              vProb[s1][s2][indx][snp] = alphaTilde[s1][indx][snp] * Tmat_ss(s1, s2, r_c[snp], r_c[snp+nSnps_c-1]) *
                Qentry(pOPGP[(snp+1)*noFam_c + fam], pAA[indx][snp+1], pAB[indx][snp+1], pBB[indx][snp+1], s2+1)  * 
                betaTilde[s2][indx][snp+1];
            }
          }
        }
        snp = nSnps_c - 1;
        for(s1 = 0; s1 < 4; s1++){
          uProb[s1][indx][snp] = (alphaTilde[s1][indx][snp] * betaTilde[s1][indx][snp])/exp(-log_w[indx][snp]); 
        }
      }  
    }

    //////// M-step:
    // The recombination fractions
    if(sexSpec_c){
      //Rprintf("Sex-Specific rf's\n");
      #pragma omp parallel for private(sum, fam, ind, s1, s2) num_threads(nThreads_c)
      for(snp = 0; snp < nSnps_c-1; snp++){
        // Paternal
        if(pss_rf[snp] == 1){
          sum = 0;
          for(fam = 0; fam < noFam_c; fam++){
            for(ind = 0; ind < nInd_c[fam]; ind++){
              int indx = ind + indSum[fam];
              for(s1 = 0; s1 < 2; s1++){
                for(s2 = 2; s2 < 4; s2++){
                  sum = sum + vProb[s1][s2][indx][snp] + vProb[s1+2][s2-2][indx][snp];
                }
              }
            }
          }
          r_c[snp] = 1.0/nTotal * sum;
        }
        // Maternal
        if(pss_rf[snp+nSnps_c-1] == 1){
          sum = 0;
          for(fam = 0; fam < noFam_c; fam++){
            for(ind = 0; ind < nInd_c[fam]; ind++){
              int indx = ind + indSum[fam];
              for(s1 = 0; s1 < 2; s1++){
                for(s2 = 0; s2 < 2; s2++){
                  sum = sum + vProb[2*s1][1+2*s2][indx][snp] + vProb[1+2*s1][2*s2][indx][snp];
                }
              }
            }
          }
          r_c[snp + nSnps_c-1] = 1.0/nTotal * sum;
        }
      }
    }
    else{ // non sex-specific
      //Rprintf("non Sex-Specific rf's\n");
      #pragma omp parallel for private(sum, fam, ind, s1, s2) num_threads(nThreads_c)
      for(snp = 0; snp < nSnps_c-1; snp++){
        sum = 0;
        for(fam = 0; fam < noFam_c; fam++){
          for(ind = 0; ind < nInd_c[fam]; ind++){
            int indx = ind + indSum[fam];
            for(s1 = 0; s1 < 4; s1++){
              for(s2 = 0; s2 < 4; s2++){
                sum = sum + vProb[s1][s2][indx][snp] * Tcount(s1,s2);
                //Rprintf("sum: %.36Le and snp %i and s1 %i and s2 %i\n", vProb[s1][s2][indx][snp] * Tcount(s1,s2), snp, s1, s2);
              }
            }
          }
        }
        //Rprintf("sum: %.36Le and snp %i\n", sum, snp);
        //Rprintf("r_c: %.36Le and snp %i\n", 1.0/(2.0*nTotal) * sum, snp);
        r_c[snp] = 1.0/(2.0*nTotal) * sum;
        r_c[snp + nSnps_c-1] = 1.0/(2.0*nTotal) * sum;
      }
    }
    // Error parameter:
    if(seqError_c){
      sumA = 0;
      sumB = 0;
      #pragma omp parallel for reduction(+:sumA,sumB) private(fam, ind, a, b, s1) num_threads(nThreads_c)
      for(snp = 0; snp < nSnps_c; snp++){
        for(fam = 0; fam < noFam_c; fam++){
          for(ind = 0; ind < nInd_c[fam]; ind++){
            int indx = ind + indSum[fam];
            a = pref[indx + nTotal*snp];
            b = palt[indx + nTotal*snp];
            //Rprintf("a :%.8f at snp %i and ind %i\n", a, snp, indx);
            //Rprintf("b :%.8f at snp %i and ind %i\n", b, snp, indx);
            for(s1 = 0; s1 < 4; s1++){
              //Rprintf("Iaa :%i at s1 %i and snp %i and ind %i\n", Iaa[s1][indx][snp], s1,snp, indx);
              //Rprintf("Ibb :%i at s1 %i and snp %i and ind %i\n", Ibb[s1][indx][snp], s1,snp, indx);
              sumA = sumA + uProb[s1][indx][snp]*(b*Iaa[s1][indx][snp] +  a*Ibb[s1][indx][snp]);
              sumB = sumB + uProb[s1][indx][snp]*(a*Iaa[s1][indx][snp] +  b*Ibb[s1][indx][snp]);

            }
          }
        }
      }
      ep_c = sumA/(sumA + sumB);
    }
    //Rprintf("llvalue :%.8f at iter %i\n", llval, iter);
    //Rprintf("prellval :%.8f at iter %i\n", prellval, iter);
    //Rprintf("diff in lik :%.8f at iter %i\n", (llval - prellval) , iter);
    //Rprintf("llvalue :%.22f at iter %i\n", llval, iter);
    //Rprintf("ep :%.22f at iter %i\n", ep_c, iter);
    //for(snp = 0; snp < nSnps_c - 1; snp++){
    //    Rprintf("r_c :%.22f at iter %i\n", r_c[snp], iter);
    //}
  }

  // Set up the R output object.
  for(snp = 0; snp < nSnps_c - 1; snp++){
    for(parent = 0; parent < 2; parent++){
      prout[snp+parent*(nSnps_c-1)] = r_c[snp+parent*(nSnps_c-1)];
    }
  }
  pepout[0] = ep_c;
  pllout[0] = llval;
  SET_VECTOR_ELT(pout, 0, rout);
  SET_VECTOR_ELT(pout, 1, epout);
  SET_VECTOR_ELT(pout, 2, llout);
  // Return the parameter estimates of log-likelihood value
  UNPROTECT(4);
  return pout;
}



// Function for extracting entries of the emission probability matrix
// when the OPGPs are known
int Iindx_up(int config, int elem){
  switch(config){
  case 1:
    if(elem == 1)
      return 2;
    else if ((elem == 2)|(elem == 3))  
      return 0;
    else if (elem == 4)
      return 1;
  case 2:
    if ((elem == 1)|(elem == 2))
      return 0;
    else if ((elem == 3)|(elem == 4))
      return 1;
  case 3:
    if ((elem == 1)|(elem == 2))
      return 2;
    else if ((elem == 3)|(elem == 4))
      return 0;
  case 4:
    if ((elem == 1)|(elem == 3))
      return 0;
    else if ((elem == 2)|(elem == 4))
      return 1;
  case 5:
    if ((elem == 1)|(elem == 3))
      return 2;
    else if ((elem == 2)|(elem == 4))
      return 0;
  } // end of Switch
  return -1;
}



SEXP EM_HMM_UP(SEXP r, SEXP ep, SEXP ref, SEXP alt, SEXP config, SEXP noFam, SEXP nInd, SEXP nSnps,
               SEXP seqError, SEXP para, SEXP ss_rf, SEXP nThreads){
  // Initialize variables
  int s1, s2, fam, ind, snp, g, iter, nIter, indx, parent, noFam_c, nSnps_c, seqError_c, nThreads_c;
  double sum, sumA, sumB, a, b, delta;

  // set up number of threads
  nThreads_c = asInteger(nThreads);
  if (nThreads_c <= 0) {
    // if nThreads is set to zero then use everything
    nThreads_c = omp_get_max_threads();
  }

  // Copy values of R input into C objects
  nSnps_c = INTEGER(nSnps)[0];
  double r_c[(nSnps_c-1)*2], ep_c;
  for(snp = 0; snp < nSnps_c - 1; snp++){
    for(parent = 0; parent < 2; parent++){
      r_c[snp+parent*(nSnps_c-1)] = REAL(r)[snp+parent*(nSnps_c-1)];
    }
  }
  ep_c = REAL(ep)[0];
  seqError_c = INTEGER(seqError)[0];
  noFam_c = INTEGER(noFam)[0];
  int nInd_c[noFam_c], nTotal, indSum[noFam_c];
  nTotal = 0;
  for(fam = 0; fam < noFam_c; fam++){
    nInd_c[fam] = INTEGER(nInd)[fam];
    indSum[fam] = nTotal;
    nTotal = nTotal + nInd_c[fam];
  }
  nIter = REAL(para)[0];
  delta = REAL(para)[1];
  // Initialize some more variables
  double alphaTilde[4][nTotal][nSnps_c], alphaDot[4];
  double betaTilde[4][nTotal][nSnps_c], betaDot[4];
  double w_new, log_w[nTotal][nSnps_c];
  double uProb[4][nTotal][nSnps_c], vProb[4][4][nTotal][nSnps_c-1];
  // Define the pointers to the other input R variables
  int *pconfig, *pref, *palt;
  pconfig = INTEGER(config);
  pref = INTEGER(ref);
  palt = INTEGER(alt);
  int *pss_rf;
  pss_rf = INTEGER(ss_rf);
  // if sex-specific rf
  for(snp = 0; snp < nSnps_c-1; snp++){
    if(pss_rf[snp]==0){
      r_c[snp] = 0;
    }
    if(pss_rf[snp+nSnps_c-1]==0){
      r_c[snp + nSnps_c-1] = 0;
    }
  }
  
  // Probability matrics
  double pAA[nTotal][nSnps_c];
  double pAB[nTotal][nSnps_c];
  double pBB[nTotal][nSnps_c];
  double bin_coef[nTotal][nSnps_c];
  int Iaa[4][nTotal][nSnps_c], Ibb[4][nTotal][nSnps_c];
  //double (*pbin_coef)[nSnps_c] = &bin_coef;
  // Compute the probs for the heterozygous calls
  for(fam = 0; fam < noFam_c; fam++){
    for(ind = 0; ind < nInd_c[fam]; ind++){
      indx = ind + indSum[fam];
      for(snp = 0; snp < nSnps_c; snp++){
        bin_coef[indx][snp] = 1; //binomial(pref[indx + nTotal*snp] + palt[indx + nTotal*snp], pref[indx + nTotal*snp]);
        pAB[indx][snp] = bin_coef[indx][snp] * powl(0.5,pref[indx + nTotal*snp] + palt[indx + nTotal*snp]);
        for(s1 = 0; s1 < 4; s1++){
          g = Iindx_up(pconfig[snp*noFam_c + fam], s1 + 1);
          if( g == 1 ){
            Iaa[s1][indx][snp] = 1;
            Ibb[s1][indx][snp] = 0;
          }
          else if( g == 2 ){
            Ibb[s1][indx][snp] = 1;
            Iaa[s1][indx][snp] = 0;
          }
          else{
            Iaa[s1][indx][snp] = 0;
            Ibb[s1][indx][snp] = 0;
          }
        }
      }
    }
  }
  // Define the output variable
  double *prout, *pepout, *pllout;
  SEXP rout = PROTECT(allocVector(REALSXP, 2*(nSnps_c-1)));
  SEXP epout = PROTECT(allocVector(REALSXP, 1));
  SEXP llout = PROTECT(allocVector(REALSXP, 1));
  prout = REAL(rout);
  pepout = REAL(epout);
  pllout = REAL(llout);
  SEXP pout = PROTECT(allocVector(VECSXP, 3));
  double llval = 0, prellval = 0;

  /////// Start algorithm
  iter = 0;
  while( (iter < 2) || ((iter < nIter) & ((llval - prellval) > delta))){
    iter = iter + 1;
    prellval = llval;
    llval = 0;
    
    // unpdate the probabilites for pAA and pBB given the parameter values and data.
    computeProb(nTotal, nSnps_c, pAA, pBB, bin_coef, ep_c, pref, palt, nThreads_c);
    
    // Compute the forward and backward probabilities for each individual
    for(fam = 0; fam < noFam_c; fam++){
      #pragma omp parallel for reduction(+:llval) num_threads(nThreads_c) \
                               private(indx, sum, s1, s2, alphaDot, snp, w_new, betaDot)
      for(ind = 0; ind < nInd_c[fam]; ind++){
        indx = ind + indSum[fam];
        //Rprintf("indSum %i\n", indSum[fam]);
        //Rprintf("indx %i\n", indx);
        //Rprintf("bin_coef %f at snp %i: ind %i:\n",bin_coef[indx][0], 0, indx);
        //Rprintf("pAA %f at snp %i: ind %i:\n",pAA[indx][0], 0, indx);
        //Rprintf("pAB %f at snp %i: ind %i:\n",pAB[indx][0], 0, indx);
        //Rprintf("pBB %f at snp %i: ind %i:\n",pBB[indx][0], 0, indx);
        /////////////////////////////////////////
        // Compute forward probabilities at snp 1
        sum = 0;
        //Rprintf("config :%.6f at snp %i in ind %i\n", pconfig[fam], 0, ind);
        for(s1 = 0; s1 < 4; s1++){
          //Rprintf("Q value :%.6f at snp %i in ind %i\n", Qentry_up(pconfig[fam], pAA[indx][0], pAB[indx][0], pBB[indx][0], s1+1), 0, ind);
          alphaDot[s1] = 0.25 * Qentry_up(pconfig[fam], pAA[indx][0], pAB[indx][0], pBB[indx][0], s1+1);
          sum = sum + alphaDot[s1];
          //Rprintf("pAA %.16f, pAB %.16f pBB %.16f at snp %i in ind %i\n", pAA[indx][0], pAB[indx][0], pBB[indx][0], 0, ind);
          //Rprintf("Q value :%.22f at snp %i in ind %i\n", Qentry_up(pconfig[fam], pAA[indx][0], pAB[indx][0], pBB[indx][0], s1+1), 0, ind);
        }
        // Scale forward probabilities
        for(s1 = 0; s1 < 4; s1++){
          alphaTilde[s1][indx][0] = alphaDot[s1]/sum;
          //Rprintf("alphaTilde :%.6f at state %i, snp %i and ind %i\n", alphaTilde[s1][indx][0], s1, 0, ind);
        }
        //Rprintf("New weight :%.6f at snp %i\n", sum, 1);
        
        // add contribution to likelihood
        log_w[indx][0] = log(sum);
        llval = llval + log(sum);
        //w_logcumsum = log(sum);
        //llval = llval + w_logcumsum;
        // llval = llval - log(sum);
        //Rprintf("llvalue :%.6f at snp %i\n", llval, 0);
        
        // iterate over the remaining SNPs
        for(snp = 1; snp < nSnps_c; snp++){
          //Rprintf("snp %i\n", snp );
          //Rprintf("nSnps %i\n", nSnps_c );
          // compute the next forward probabilities for snp j
          //Rprintf("pAA %f at snp %i: ind %i:\n",pAA[indx][snp], snp, ind);
          //Rprintf("pAB %f at snp %i: ind %i:\n",pAB[indx][snp], snp, ind);
          //Rprintf("pBB %f at snp %i: ind %i:\n",pBB[indx][snp], snp, ind);
          for(s2 = 0; s2 < 4; s2++){
            sum = 0;
            for(s1 = 0; s1 < 4; s1++){
              sum = sum + Tmat_ss(s1, s2, r_c[snp-1], r_c[snp-1+nSnps_c-1]) * alphaTilde[s1][indx][snp-1];
              //Rprintf("Tmat %.16f", Tmat_ss(s1, s2, r_c[snp-1], r_c[snp-1+nSnps_c-1]));
            }
            //Rprintf("sum :%.16f at state %i, snp %i and ind %i\n", sum, s2, snp, ind);
            //Rprintf("pAA %.16f, pAB %.16f pBB %.16f at snp %i in ind %i\n", pAA[indx][snp], pAB[indx][snp], pBB[indx][snp], snp, indx);
            //Rprintf("Q value :%.16f at snp %i in ind %i\n", Qentry_up(pconfig[snp*noFam_c + fam], pAA[indx][snp], pAB[indx][snp], pBB[indx][snp], s2+1), snp, ind);
            alphaDot[s2] = Qentry_up(pconfig[snp*noFam_c + fam], pAA[indx][snp], pAB[indx][snp], pBB[indx][snp], s2+1) * sum;
            //Rprintf("alphaDot :%.22f at state %i, snp %i and ind %i\n", alphaDot[s2], s2, snp, ind);
          }
          
          // Compute the weight for snp \ell
          w_new = 0;
          for(s2 = 0; s2 < 4; s2++){
            w_new = w_new + alphaDot[s2];
          }
          // Add contribution to the likelihood
          //llval = llval + log(w_new) + w_logcumsum;
          //w_logcumsum = w_logcumsum + log(w_new);
          log_w[indx][snp] = log(w_new);
          // Scale the forward probability vector
          for(s2 = 0; s2 < 4; s2++){
            alphaTilde[s2][indx][snp] = alphaDot[s2]/w_new;
          }
          // Add contribution to the likelihood
          llval = llval + log(w_new);
          //Rprintf("llvalue :%.8f at ind %i\n", llval, indx);
        }
        //////////////////////
        // Compute the backward probabilities
        for(s1 = 0; s1 < 4; s1++){
          betaTilde[s1][indx][nSnps_c-1] = 1/exp(log_w[indx][nSnps_c-1]);
          //sum = log(betaTilde[s1][indx][nSnps_c-1]) - (log(1) -  log_w[indx][nSnps_c-1]);
          //Rprintf("log(beta) :%.22f at state %i, snp %i and ind %i\n", exp(sum), s1, nSnps_c-1, ind);
        }
        // iterate over the remaining SNPs
        sum = 0;
        for(snp = nSnps_c-2; snp > -1; snp--){
          for(s1 = 0; s1 < 4; s1++){
            sum = 0;
            for(s2 = 0; s2 < 4; s2++){
              sum = sum + Qentry_up(pconfig[(snp+1)*noFam_c + fam], pAA[indx][snp+1], pAB[indx][snp+1], pBB[indx][snp+1], s2+1) * 
                Tmat_ss(s1, s2, r_c[snp], r_c[snp+nSnps_c-1]) * betaTilde[s2][indx][snp+1];
            }
            betaDot[s1] = sum;
          }
          // Scale the backward probability vector
          for(s1 = 0; s1 < 4; s1++){
            betaTilde[s1][indx][snp] = betaDot[s1]/exp(log_w[indx][snp]);
            //Rprintf("betaTilde :%.6f at state %i, snp %i and ind %i\n", betaTilde[s1][indx][snp], s1, snp, ind);
          }
        }
        ///////// E-step:
        for(snp = 0; snp < nSnps_c - 1; snp++){
          for(s1 = 0; s1 < 4; s1++){
            uProb[s1][indx][snp] = (alphaTilde[s1][indx][snp] * betaTilde[s1][indx][snp])/exp(-log_w[indx][snp]);
            for(s2 = 0; s2 < 4; s2++){
              vProb[s1][s2][indx][snp] = alphaTilde[s1][indx][snp] * Tmat_ss(s1, s2, r_c[snp], r_c[snp+nSnps_c-1]) *
                Qentry_up(pconfig[(snp+1)*noFam_c + fam], pAA[indx][snp+1], pAB[indx][snp+1], pBB[indx][snp+1], s2+1)  * 
                betaTilde[s2][indx][snp+1];
            }
          }
        }
        snp = nSnps_c - 1;
        for(s1 = 0; s1 < 4; s1++){
          uProb[s1][indx][snp] = (alphaTilde[s1][indx][snp] * betaTilde[s1][indx][snp])/exp(-log_w[indx][snp]); 
        }
      }  
    }
    // The recombination fractions
    #pragma omp parallel for private(sum, fam, ind, indx, s1, s2) num_threads(nThreads_c)
    for(snp = 0; snp < nSnps_c-1; snp++){
      // Paternal
      if(pss_rf[snp] == 1){
        sum = 0;
        for(fam = 0; fam < noFam_c; fam++){
          for(ind = 0; ind < nInd_c[fam]; ind++){
            indx = ind + indSum[fam];
            for(s1 = 0; s1 < 2; s1++){
              for(s2 = 2; s2 < 4; s2++){
                sum = sum + vProb[s1][s2][indx][snp] + vProb[s1+2][s2-2][indx][snp];
              }
            }
          }
        }
        r_c[snp] = 1.0/nTotal * sum;
      }
      // Maternal
      if(pss_rf[snp+nSnps_c-1] == 1){
        sum = 0;
        for(fam = 0; fam < noFam_c; fam++){
          for(ind = 0; ind < nInd_c[fam]; ind++){
            indx = ind + indSum[fam];
            for(s1 = 0; s1 < 2; s1++){
              for(s2 = 0; s2 < 2; s2++){
                sum = sum + vProb[2*s1][1+2*s2][indx][snp] + vProb[1+2*s1][2*s2][indx][snp];
              }
            }
          }
        }
        r_c[snp + nSnps_c-1] = 1.0/nTotal * sum;
      }
    }
    // Error parameter:
    if(seqError_c){
      sumA = 0;
      sumB = 0;
      #pragma omp parallel for reduction(+:sumA,sumB) num_threads(nThreads_c) \
                               private(fam, ind, indx, a, b, s1)
      for(snp = 0; snp < nSnps_c; snp++){
        for(fam = 0; fam < noFam_c; fam++){
          for(ind = 0; ind < nInd_c[fam]; ind++){
            indx = ind + indSum[fam];
            a = pref[indx + nTotal*snp];
            b = palt[indx + nTotal*snp];
            for(s1 = 0; s1 < 4; s1++){
              sumA = sumA + uProb[s1][indx][snp]*(b*Iaa[s1][indx][snp] +  a*Ibb[s1][indx][snp]);
              sumB = sumB + uProb[s1][indx][snp]*(a*Iaa[s1][indx][snp] +  b*Ibb[s1][indx][snp]);
            }
          }
        }
      }
      ep_c = sumA/(sumA + sumB);
    }
    //Rprintf("llvalue :%.8f at iter %i\n", llval, iter);
    //Rprintf("prellval :%.8f at iter %i\n", prellval, iter);
    //Rprintf("diff in lik :%.8f at iter %i\n", (llval - prellval) , iter);
    //Rprintf("ep :%.22f at iter %i\n", ep_c, iter);
    //for(snp = 0; snp < nSnps_c - 1; snp++){
    //  Rprintf("r_c :%.22f at iter %i\n", r_c[snp], iter);
    //}
  }

  
  // Set up the R output object.
  for(snp = 0; snp < nSnps_c - 1; snp++){
    for(parent = 0; parent < 2; parent++){
      prout[snp+parent*(nSnps_c-1)] = r_c[snp+parent*(nSnps_c-1)];
    }
  }
  pepout[0] = ep_c;
  pllout[0] = llval;
  SET_VECTOR_ELT(pout, 0, rout);
  SET_VECTOR_ELT(pout, 1, epout);
  SET_VECTOR_ELT(pout, 2, llout);
  // Return the parameter estimates of log-likelihood value
  UNPROTECT(4);
  return pout;
}



