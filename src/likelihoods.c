#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

// Function for extracting entries of the emission probability matrix
// when the OPGPs are known
double Qentry(int OPGP,double Kaa,double Kab, double Kbb,int elem){
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
// when the OPGP are considered the baseline (and so phase is unknown and the r.f's are sex-specific)
double Qentry_up(int config,double Kaa,double Kab, double Kbb,int elem){
  switch(config){
  case 1:
    if(elem == 1)
      return Kbb;
    else if ((elem == 2)|(elem == 3))  
      return Kab;
    else if (elem == 4)
      return Kaa;
  case 2:
    if ((elem == 1)|(elem == 2))
      return Kab;
    else if ((elem == 3)|(elem == 4))
      return Kaa;
  case 3:
    if ((elem == 1)|(elem == 2))
      return Kbb;
    else if ((elem == 3)|(elem == 4))
      return Kab;
  case 4:
    if ((elem == 1)|(elem == 3))
      return Kab;
    else if ((elem == 2)|(elem == 4))
      return Kaa;
  case 5:
    if ((elem == 1)|(elem == 3))
      return Kbb;
    else if ((elem == 2)|(elem == 4))
      return Kab;
  } // end of Switch
  return -1;
}


// Function for returning a specified enetry of the transition matrix for a given recombination fraction value
double Tmat(int s1, int s2, double rval){
  int sSum = s1 + s2*4;
  if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15))
    return (1-rval)*(1-rval);
  else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12))
    return rval*rval;
  else
    return (1-rval)*rval;
}

// Function for returning a specified enetry of the transition matrix for a given recombination fraction value
// when the r.f.'s are sex-specific 
double Tmat_ss(int s1, int s2, double r_f, double r_m){
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
SEXP ll_fs_scaled_err_c(SEXP r, SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP OPGP, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int s1, s2, ind, snp, nInd_c, nSnps_c;
  double *pll, *pr, *pKaa, *pKab, *pKbb, *pOPGP;
  double alphaTilde[4], alphaDot[4], sum, w_logcumsum, w_new;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = REAL(OPGP);
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
    w_logcumsum = log(sum);
    llval = llval + w_logcumsum;

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
      llval = llval + log(w_new) + w_logcumsum;
      w_logcumsum = w_logcumsum + log(w_new);
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
  int s1, s2, ind, snp, nInd_c, nSnps_c;
  double *pll, *pr, *pKaa, *pKab, *pKbb, *pOPGP;
  double alphaTilde[4], alphaDot[4], sum, w_logcumsum, w_new;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = REAL(OPGP);
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
    w_logcumsum = log(sum);
    llval = llval + w_logcumsum;

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
      llval = llval + log(w_new) + w_logcumsum;
      w_logcumsum = w_logcumsum + log(w_new);
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
  int s1, s2, ind, snp, nInd_c, nSnps_c;
  double *pll, *pr, *pKaa, *pKab, *pKbb, *pconfig;
  double alphaTilde[4], alphaDot[4], sum, w_logcumsum, w_new;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pconfig = REAL(config);
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
    w_logcumsum = log(sum);
    llval = llval + w_logcumsum;
    
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
      llval = llval + log(w_new) + w_logcumsum;
      w_logcumsum = w_logcumsum + log(w_new);
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
