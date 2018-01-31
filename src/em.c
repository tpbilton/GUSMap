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

// Function for returning a specified enetry of the transition matrix for a given recombination fraction value
// when the r.f.'s are sex-specific 
long double Tmat_ss(int s1, int s2, long double r_f, long double r_m){
  int sSum = s1 + s2*4;
  if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15))
    return (1-r_f)*(1-r_m);
  else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12))
    return r_f*r_m;
  else if((sSum == 1)|(sSum == 4)|(sSum == 11)|(sSum == 14))
    return (1-r_f)*r_m;
  else 
    return r_f*(1-r_m);
}

int Tcount(int s1, int s2){
  int sSum = s1 + s2*4;
  if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15))
    return 0;
  else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12))
    return 2;
  else
    return 1;
}

// Function for computing binomial coefficients
// taken  from this website "https://rosettacode.org/wiki/Evaluate_binomial_coefficients#C"
static unsigned long gcd_ui(unsigned long x, unsigned long y) {
  unsigned long t;
  if (y < x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
}

unsigned long binomial(unsigned long n, unsigned long k) {
  unsigned long d, g, r = 1;
  if (k == 0) return 1;
  if (k == 1) return n;
  if (k >= n) return (k == n);
  if (k > n/2) k = n-k;
  for (d = 1; d <= k; d++) {
    if (r >= ULONG_MAX/n) {  /* Possible overflow */
unsigned long nr, dr;  /* reduced numerator / denominator */
g = gcd_ui(n, d);  nr = n/g;  dr = d/g;
g = gcd_ui(r, dr);  r = r/g;  dr = dr/g;
if (r >= ULONG_MAX/nr) return 0;  /* Unavoidable overflow */
r *= nr;
r /= dr;
n--;
    } else {
      r *= n--;
      r /= d;
    }
  }
  return r;
}


// Function for computing the emission probabilities given the true genotypes
long double computeProb(long double *ppAA, long double *ppBB, long double *pbin_coef,
                        long double epsilon, double *pdepth_Ref, double *pdepth_Alt,
                        int nInd, int nSnps){
  int snp, ind, indx;
  for(ind = 0; ind < nInd; ind++){
    for(snp = 0; snp < nSnps; snp++){
      indx = ind + nInd * snp;
      if( (pdepth_Ref[indx] + pdepth_Alt[indx]) == 0){
        *ppAA = 1;
        *ppBB = 1;
      }
      else{
        *ppAA = *pbin_coef * powl(1 - epsilon, pdepth_Ref[indx]) * powl(epsilon, pdepth_Alt[indx]);
        *ppBB = *pbin_coef * powl(epsilon, pdepth_Ref[indx]) * powl(1 - epsilon, pdepth_Alt[indx]);
      }
      ppAA++;
      ppBB++;
      pbin_coef++;
    }
  }
  return -1;
}

// Function for extracting entries of the emission probability matrix
// when the OPGPs are known
long double Qentry(int OPGP, long double Kaa, long double Kab, long double Kbb,int elem){
  switch(OPGP){
  case 1:
    if(elem == 1)
      return Kbb;
    else if ((elem == 2)|(elem == 3))  
      return Kab;
    else if (elem == 4)
      return Kaa;
  case 2:
    if(elem == 3)
      return Kbb;
    else if ((elem == 1)|(elem == 4))
      return Kab;
    else if (elem == 2)
      return Kaa;
  case 3:
    if(elem == 2) 
      return Kbb;
    else if ((elem == 1)|(elem == 4))
      return Kab;
    else if (elem == 3)
      return Kaa;
  case 4:
    if(elem == 4) 
      return Kbb;
    else if ((elem == 2)|(elem == 3))
      return Kab;
    else if (elem == 1)
      return Kaa;
  case 5:
    if ((elem == 1)|(elem == 2))
      return Kab;
    else if ((elem == 3)|(elem == 4))
      return Kaa;
  case 6:
    if ((elem == 1)|(elem == 2))
      return Kaa;
    else if ((elem == 3)|(elem == 4))
      return Kab;
  case 7:
    if ((elem == 1)|(elem == 2))
      return Kbb;
    else if ((elem == 3)|(elem == 4))
      return Kab;
  case 8:
    if ((elem == 1)|(elem == 2))
      return Kab;
    else if ((elem == 3)|(elem == 4))
      return Kbb;
  case 9:
    if ((elem == 1)|(elem == 3))
      return Kab;
    else if ((elem == 2)|(elem == 4))
      return Kaa;
  case 10:
    if ((elem == 1)|(elem == 3))
      return Kaa;
    else if ((elem == 2)|(elem == 4))
      return Kab;
  case 11:
    if ((elem == 1)|(elem == 3))
      return Kbb;
    else if ((elem == 2)|(elem == 4))
      return Kab;
  case 12:
    if ((elem == 1)|(elem == 3))
      return Kab;
    else if ((elem == 2)|(elem == 4))
      return Kbb;
  case 13:
    return Kaa;
  case 14:
    return Kab;
  case 15:
    return Kab;
  case 16:
    return Kbb;
  } // end of Switch
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




SEXP EM_HMM(SEXP r, SEXP ep, SEXP depth_Ref, SEXP depth_Alt, SEXP OPGP, SEXP noFam, SEXP nInd, SEXP nSnps,
            SEXP sexSpec, SEXP seqError, SEXP para){
  // Initialize variables
  int s1, s2, fam, ind, snp, g, iter, nIter, indx, parent, noFam_c, nSnps_c, sexSpec_c, seqError_c;
  long double sum, sumA, sumB, a, b, delta;
  // Copy values of R input into C objects
  nSnps_c = INTEGER(nSnps)[0];
  long double r_c[(nSnps_c-1)*2], ep_c;
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
    nTotal = nTotal + nInd_c[fam];
    indSum[fam] = nTotal - nInd_c[fam];
  }
  nIter = REAL(para)[0];
  delta = REAL(para)[1];
  // Initialize some more variables
  long double alphaTilde[4][nTotal][nSnps_c], alphaDot[4];
  long double betaTilde[4][nTotal][nSnps_c], betaDot[4];
  long double w_logcumsum, w_new, log_w[nTotal][nSnps_c];
  long double uProb[4][nTotal][nSnps_c], vProb[4][4][nTotal][nSnps_c-1];
  // Define the pointers to the other input R variables
  double *pOPGP, *pdepth_Ref, *pdepth_Alt;
  pOPGP = REAL(OPGP);
  pdepth_Ref = REAL(depth_Ref);
  pdepth_Alt = REAL(depth_Alt);
  // Probability matrics
  long double pAA[nTotal][nSnps_c];
  long double pAB[nTotal][nSnps_c];
  long double pBB[nTotal][nSnps_c];
  long double bin_coef[nTotal][nSnps_c];
  long double *ppAA, *ppBB, *pbin_coef;
  int Iaa[4][nTotal][nSnps_c], Ibb[4][nTotal][nSnps_c];
  ppAA = &pAA[0][0];
  ppBB = &pBB[0][0];
  pbin_coef = &bin_coef[0][0];
  //double (*pbin_coef)[nSnps_c] = &bin_coef;
  // Compute the probs for the heterozygous calls
  for(fam = 0; fam < noFam_c; fam++){
    for(ind = 0; ind < nInd_c[fam]; ind++){
      indx = ind + indSum[fam];
      for(snp = 0; snp < nSnps_c; snp++){
        bin_coef[indx][snp] = binomial(pdepth_Ref[indx + nTotal*snp] + pdepth_Alt[indx + nTotal*snp], pdepth_Ref[indx + nTotal*snp]);
        pAB[indx][snp] = bin_coef[indx][snp] * powl(0.5,pdepth_Ref[indx + nTotal*snp] + pdepth_Alt[indx + nTotal*snp]);
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
  long double llval = 0, prellval = 0;
  Rprintf("delta :%.22f\n", delta);
  
  /////// Start algorithm
  iter = 0;
  while( (iter < 2) || ((iter < nIter) & ((llval - prellval) > delta))){
    iter = iter + 1;
    prellval = llval;
    llval = 0;
    
    // unpdate the probabilites for pAA and pBB given the parameter values and data.
    computeProb(ppAA, ppBB, pbin_coef, ep_c, pdepth_Ref, pdepth_Alt, nTotal, nSnps_c);
  
    // Compute the forward and backward probabilities for each individual
    for(fam = 0; fam < noFam_c; fam++){
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
        //Rprintf("OPGP :%.6f at snp %i in ind %i\n", pOPGP[fam], 0, ind);
        for(s1 = 0; s1 < 4; s1++){
          //Rprintf("Q value :%.6f at snp %i in ind %i\n", Qentry(pOPGP[fam], pAA[indx][0], pAB[indx][0], pBB[indx][0], s1+1), 0, ind);
          alphaDot[s1] = 0.25 * Qentry(pOPGP[fam], pAA[indx][0], pAB[indx][0], pBB[indx][0], s1+1);
          sum = sum + alphaDot[s1];
        }
        // Scale forward probabilities
        for(s1 = 0; s1 < 4; s1++){
          alphaTilde[s1][indx][0] = alphaDot[s1]; ///sum;
          //Rprintf("alphaTilde :%.6f at state %i, snp %i and ind %i\n", alphaTilde[s1][indx][0], s1, 0, ind);
        }
        //Rprintf("New weight :%.6f at snp %i\n", sum, 1);
      
        // add contribution to likelihood
        log_w[indx][0] = log(1); //log(sum);
        w_logcumsum = log(sum);
        llval = llval + w_logcumsum;
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
            }
            //Rprintf("sum :%.6f at state %i, snp %i and ind %i\n", sum, s2, snp, ind);
            //Rprintf("Q value :%.6f at snp %i in ind %i\n", Qentry(pOPGP[snp*noFam_c + fam], pAA[indx][snp], pAB[indx][snp], pBB[indx][snp], s2+1), snp, ind);
            alphaDot[s2] = Qentry(pOPGP[snp*noFam_c + fam], pAA[indx][snp], pAB[indx][snp], pBB[indx][snp], s2+1) * sum;
            //Rprintf("alphaDot :%.22f at state %i, snp %i and ind %i\n", alphaDot[s2], s2, snp, ind);
          }
          
          // Compute the weight for snp \ell
          w_new = 0;
          for(s2 = 0; s2 < 4; s2++){
            w_new = w_new + alphaDot[s2];
          }
          // Add contribution to the likelihood
          llval = llval + log(w_new) + w_logcumsum;
          w_logcumsum = w_logcumsum + log(w_new);
          log_w[indx][snp] = log(1); //log(w_new);
          // Scale the forward probability vector
          for(s2 = 0; s2 < 4; s2++){
            alphaTilde[s2][indx][snp] = alphaDot[s2];// /w_new;
          }
        }
        //////////////////////
        // Compute the backward probabilities
        for(s1 = 0; s1 < 4; s1++){
          betaTilde[s1][indx][nSnps_c-1] = 1;//exp(-log_w[indx][nSnps_c-1]);
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
            betaTilde[s1][indx][snp] = betaDot[s1]; ///exp(log_w[indx][snp]);
            //Rprintf("betaTilde :%.6f at state %i, snp %i and ind %i\n", betaTilde[s1][indx][snp], s1, snp, ind);
          }
        }
      }  
    }
/*    fam = 1;
    ind = 1;
    indx = ind + indSum[fam];
    sumA = 0;
    for(snp = 0; snp < nSnps_c; snp++){
      sumA = sumA + log_w[indx][snp];
      Rprintf("sumA :%.22f at snp %i and ind %i\n", log_w[indx][snp], snp, ind);
      for(s1 = 0; s1 < 4; s1++){
        Rprintf("alphaTilde :%.22f at state %i, snp %i and ind %i\n", alphaTilde[s1][indx][snp], s1, snp, ind);
      }
    }
    sumB = 0;
    for(snp = nSnps_c-1; snp > -1; snp--){
      sumB = sumB + log_w[indx][snp];
      Rprintf("sumB :%.22f at snp %i and ind %i\n", log_w[indx][snp], snp, ind);
      for(s1 = 0; s1 < 4; s1++){
        Rprintf("betaTilde :%.22f at state %i, snp %i and ind %i\n", betaTilde[s1][indx][snp], s1, snp, ind);
      }
    }  */
    
    ///////// E-step:
    // Compute the probabilities of the hidden variables given the data and parameter values
    for(fam = 0; fam < noFam_c; fam++){
      for(ind = 0; ind < nInd_c[fam]; ind++){
        indx = ind + indSum[fam];
        for(snp = 0; snp < nSnps_c - 1; snp++){
          for(s1 = 0; s1 < 4; s1++){
            uProb[s1][indx][snp] = (alphaTilde[s1][indx][snp] * betaTilde[s1][indx][snp])/exp(llval); //   exp(-log_w[indx][snp]);
            for(s2 = 0; s2 < 4; s2++){
              vProb[s1][s2][indx][snp] = alphaTilde[s1][indx][snp] * Tmat_ss(s1, s2, r_c[snp], r_c[snp+nSnps_c-1]) *
                Qentry(pOPGP[(snp+1)*noFam_c + fam], pAA[indx][snp+1], pAB[indx][snp+1], pBB[indx][snp+1], s2+1)  * 
                betaTilde[s2][indx][snp+1] / exp(llval);
            }
          }
        }
        snp = nSnps_c - 1;
        for(s1 = 0; s1 < 4; s1++){
          uProb[s1][indx][snp] = (alphaTilde[s1][indx][snp] * betaTilde[s1][indx][snp])/exp(llval); //exp(-log_w[indx][snp]); 
        }
      }
    }
    for(fam = 0; fam < noFam_c; fam++){
      for(ind = 0; ind < nInd_c[fam]; ind++){
        indx = ind + indSum[fam];
        for(snp = 0; snp < nSnps_c-1; snp++){
          for(s1 = 0; s1 < 4; s1++){
            if(uProb[s1][indx][snp]>1){
              //Rprintf("uProbs: %.36Le at s1 %i and snp %i and ind %i\n", uProb[s1][indx][snp], s1, snp, ind);
            }
            for(s2 = 0; s2 < 4; s2++){
              if(vProb[s1][s2][indx][snp]>1){
                //Rprintf("vProbs: %.36Le at s1 %i and s2 %i and snp %i and ind %i\n", vProb[s1][s2][indx][snp], s1, snp, ind);
              }
              //Rprintf("vProbs: %.22f at s1 %i and s2 %i and snp %i and ind %i\n", vProb[s1][s2][indx][snp], s1, s2, snp, ind);
            }
          }
        }
      }
    }
    //////// M-step:
    // The recombination fractions
    if(sexSpec_c){
      //Rprintf("Sex-Specific rf's\n");
      for(snp = 0; snp < nSnps_c-1; snp++){
        // Paternal
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
        // Maternal
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
    else{ // non sex-specific
      //Rprintf("non Sex-Specific rf's\n");
      for(snp = 0; snp < nSnps_c-1; snp++){
        sum = 0;
        for(fam = 0; fam < noFam_c; fam++){
          for(ind = 0; ind < nInd_c[fam]; ind++){
            indx = ind + indSum[fam];
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
      for(snp = 0; snp < nSnps_c; snp++){
        for(fam = 0; fam < noFam_c; fam++){
          for(ind = 0; ind < nInd_c[fam]; ind++){
            indx = ind + indSum[fam];
            for(s1 = 0; s1 < 4; s1++){
              a = pdepth_Ref[indx + nTotal*snp];
              b = pdepth_Alt[indx + nTotal*snp];
              sumA = sumA + uProb[s1][indx][snp]*(b*Iaa[s1][indx][snp] +  a*Ibb[s1][indx][snp]);
              sumB = sumB + uProb[s1][indx][snp]*(a*Iaa[s1][indx][snp] +  b*Ibb[s1][indx][snp]);
            }
          }
        }
      }
      ep_c = sumA/(sumA + sumB);
    }
    else{
      ep_c = 0;
    }
    Rprintf("llvalue :%.36Le at iter %i\n", llval, iter);
    Rprintf("prellval :%.20Le at iter %i\n", prellval, iter);
    Rprintf("diff in lik :%.20Le at iter %i\n", (llval - prellval) , iter);
    //Rprintf("llvalue :%.22f at iter %i\n", llval, iter);
    //Rprintf("ep :%.36Le at iter %i\n", ep_c, iter);
    //for(snp = 0; snp < nSnps_c - 1; snp++){
    //    Rprintf("r_c :%.36Le at iter %i\n", r_c[snp], iter);
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


