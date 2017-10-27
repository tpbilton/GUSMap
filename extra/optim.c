#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <limits.h>


//void LevMarqAlg()

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
                        long double epsilon, int *pdepth_Ref, int *pdepth_Alt,
                        int nInd, int nSnps){
  int snp, ind, indx;
  for(snp = 0; snp < nSnps; snp++){
    for(ind = 0; ind < nInd; ind++){
      indx = ind + nInd * snp;
      if( (pdepth_Ref[indx] + pdepth_Alt[indx]) == 0){
        *ppAA = 1;
        *ppBB = 1;
      }
      else{
        *ppAA = pbin_coef[indx] * powl(1 - epsilon, pdepth_Ref[indx]) * powl(epsilon, pdepth_Alt[indx]);
        *ppBB = pbin_coef[indx] * powl(epsilon, pdepth_Ref[indx]) * powl(1 - epsilon, pdepth_Alt[indx]);
      }
      ppAA++;
      ppBB++;
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

// Function for returning a specified entry of the transition matrix for a given recombination fraction value
long double Tmat(int s1, int s2, double rval){
  int sSum = s1 + s2*4;
  if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15))
    return (1-rval)*(1-rval);
  else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12))
    return rval*rval;
  else
    return (1-rval)*rval;
}

// Function for computing the likelihood, gradient and Hessian
void ll_fs_scaled_err(long double *pll, long double *pr, long double epsilon, long double *ppAA, long double *ppAB, long double *ppBB,
                 int *pOPGP, int nInd, int nSnps){
  
  int s1, s2, ind, snp, indx;
  long double alphaTilde[4], alphaDot[4], sum, w_logcumsum, w_new, llval;
  llval = 0;
  
  // Compute the likelihood
  for(ind = 0; ind < nInd; ind++){
    // Compute forward probabilities at snp 1
    sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      Rprintf("Q value :%Le at snp %i in ind %i\n", Qentry(pOPGP[0], ppAA[ind], ppAB[ind], ppBB[ind], s1+1),0,ind);
      alphaDot[s1] = 0.25 * Qentry(pOPGP[0], ppAA[ind], ppAB[ind], ppBB[ind], s1+1);
      sum = sum + alphaDot[s1];
    }
    // Scale forward probabilities
    for(s1 = 0; s1 < 4; s1++){
      alphaTilde[s1] = alphaDot[s1]/sum;
    }
    //Rprintf("New weight :%.6f at snp %i\n", sum, 1);
    
    // Derivatives at snp 1
    
    
    // add contribution to likelihood
    w_logcumsum = log(sum);
    llval = llval + w_logcumsum;
    
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps; snp++){
      // compute the next forward probabilities for snp \ell
      for(s2 = 0; s2 < 4; s2++){
        sum = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum = sum + Tmat(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        indx = ind + nInd * snp;
        Rprintf("Q value :%Le at snp %i in ind %i\n", Qentry(pOPGP[snp], ppAA[indx], ppAB[indx], ppBB[indx], s2+1), snp, ind);
        alphaDot[s2] = Qentry(pOPGP[snp], ppAA[indx], ppAB[indx], ppBB[indx], s2+1) * sum;
      }
      // Compute the weight for snp \ell
      w_new = 0;
      for(s2 = 0; s2 < 4; s2++){
        w_new = w_new + alphaDot[s2];
      }
      // Add contribution to the likelihood
      llval = llval + log(w_new) + w_logcumsum;
      Rprintf("like value:%Le\n", llval);
      w_logcumsum = w_logcumsum + log(w_new);
      // Scale the forward probability vector
      for(s2 = 0; s2 < 4; s2++){
        alphaTilde[s2] = alphaDot[s2]/w_new;
      }
    }
  }
  *pll = -1*llval;
}
  



// Function for reading in data from R
SEXP optim_ll(SEXP r, SEXP epsilon, SEXP depth_Ref, SEXP depth_Alt, SEXP depth, SEXP OPGP, SEXP nInd, SEXP nSnps){
  
  // Variables to be copied
  int nInd_c, nSnps_c, snp, ind, indx;
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  double temp;
  long double epsilon_c, r_c[nSnps_c-1], *pr_c;
  temp = REAL(epsilon)[0];
  epsilon_c = temp;
  for(snp = 0; snp < nSnps_c - 1; snp++){
    temp = REAL(r)[snp];
    r_c[snp] = temp;
  }
  pr_c = &r_c[0];
  
  // R variables with pointers
  int *pdepth_Ref, *pdepth_Alt, *pdepth, *pOPGP;
  pdepth_Ref = INTEGER(depth_Ref);
  pdepth_Alt = INTEGER(depth_Alt);
  pdepth = INTEGER(depth);
  pOPGP = INTEGER(OPGP);
  
  // Probability variables
  long double pAA[nInd_c * nSnps_c], pAB[nInd_c * nSnps_c], pBB[nInd_c * nSnps_c], bin_coef[nInd_c * nSnps_c];
  long double *ppAA, *ppAB, *ppBB, *pbin_coef;
  ppAA = &pAA[0];
  ppAB = &pAB[0];
  ppBB = &pBB[0];
  pbin_coef = &bin_coef[0];
  
  // Define the output variable
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  double *pll;
  pll = REAL(ll);
  long double ll_c, *pll_c;
  ll_c = 0;
  pll_c = &ll_c;
  
  // Compute the values for the pAB array from the given data
  for(snp = 0; snp < nSnps_c; snp++){
    for(ind = 0; ind < nInd_c; ind++){
      indx = ind + nInd_c * snp;
      if (pdepth[indx] == 0){
        bin_coef[indx] = 1;
        pAB[indx] = 1;
      }
      else{
        bin_coef[indx] = binomial(pdepth[indx], pdepth_Ref[indx]);
        pAB[indx] = bin_coef[indx] * powl(0.5,pdepth[indx]);
      }
    }
  }
  computeProb(ppAA, ppBB, pbin_coef, epsilon_c, pdepth_Ref, pdepth_Alt, nInd_c, nSnps_c);
  
  ll_fs_scaled_err(pll_c, pr_c, epsilon_c, ppAA, ppAB, ppBB, pOPGP, nInd_c, nSnps_c);
  Rprintf("Likelihood value of function is: %Le\n", *pll_c);
  Rprintf("Likelihood value of function is: %Le\n", ll_c);
    
  for(ind = 0; ind < nInd_c; ind++){
    for(snp = 0; snp < nSnps_c; snp++){
      indx = ind + nInd_c * snp;
      Rprintf("pAA at ind %i and snp %i is: %Le\n", ind, snp, pAA[ind + nInd_c * snp]);
      Rprintf("pBB at ind %i and snp %i is: %Le\n", ind, snp, pBB[ind + nInd_c * snp]);
    }
  }
  

  
  for(snp = 0; snp < nSnps_c; snp++)
    Rprintf("Value of OPGP vector at indices %i is: %i\n", snp, pOPGP[snp]);
  
  pll[0] = ll_c;
  // Clean up and return likelihood value
  UNPROTECT(1);
  return ll;
}
  
    
  