#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

// Function for extracting entries of the emission probability matrix
// when the OPGPs are known
double Qentry(int OPGP,int g,int d,int elem, double epsilon, double delta){
  // Check first whether the genotype is missing. If so, return 1
  if (g == 4)
    return 1;
  else{ // If not missing proceed to compute value
    switch(OPGP){
    case 1:
      if(g == 1){ // AA*
        if(elem == 1) 
          return 0;
        else if ((elem == 2)|(elem == 3))
          return pow(0.5,d);
        else if (elem == 4)
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 1)|(elem == 4))
          return 0;
        else if ((elem == 2)|(elem == 3))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if (elem == 4)
          return 0;
        else if ((elem == 2)|(elem == 3))
          return pow(0.5,d);
        else if (elem == 1)
          return 1 ;
      }
    case 2:
      if(g == 1){ // AA*
        if(elem == 3) 
          return 0;
        else if ((elem == 1)|(elem == 4))
          return pow(0.5,d);
        else if (elem == 2)
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 2)|(elem == 3))
          return 0;
        else if ((elem == 1)|(elem == 4))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if (elem == 2)
          return 0;
        else if ((elem == 1)|(elem == 4))
          return pow(0.5,d);
        else if (elem == 3)
          return 1;
      }
    case 3:
      if(g == 1){ // AA*
        if(elem == 2) 
          return 0;
        else if ((elem == 1)|(elem == 4))
          return pow(0.5,d);
        else if (elem == 3)
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 2)|(elem == 3))
          return 0;
        else if ((elem == 1)|(elem == 4))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if (elem == 3)
          return 0;
        else if ((elem == 1)|(elem == 4))
          return pow(0.5,d);
        else if (elem == 2)
          return 1 ;
      }
    case 4:
      if(g == 1){ // AA*
        if(elem == 4) 
          return 0;
        else if ((elem == 2)|(elem == 3))
          return pow(0.5,d);
        else if (elem == 1)
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 1)|(elem == 4))
          return 0;
        else if ((elem == 2)|(elem == 3))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if (elem == 1)
          return 0;
        else if ((elem == 2)|(elem == 3))
          return pow(0.5,d);
        else if (elem == 4)
          return 1 ;
      }
    case 5:
      if(g == 1){ // AA*
        if ((elem == 1)|(elem == 2))
          return pow(0.5,d);
        else if ((elem == 3)|(elem == 4))
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 3)|(elem == 4))
          return 0;
        else if ((elem == 1)|(elem == 2))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if ((elem == 3)|(elem == 4))
          return 0;
        else if ((elem == 1)|(elem == 2))
          return pow(0.5,d);
      }
    case 6:
      if(g == 1){ // AA*
        if ((elem == 3)|(elem == 4))
          return pow(0.5,d);
        else if ((elem == 1)|(elem == 2))
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 1)|(elem == 2))
          return 0;
        else if ((elem == 3)|(elem == 4))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if ((elem == 1)|(elem == 2))
          return 0;
        else if ((elem == 3)|(elem == 4))
          return pow(0.5,d);
      }
    case 7:
      if(g == 1){ // AA*
        if ((elem == 1)|(elem == 3))
          return pow(0.5,d);
        else if ((elem == 2)|(elem == 4))
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 2)|(elem == 4))
          return 0;
        else if ((elem == 1)|(elem == 3))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if ((elem == 2)|(elem == 4))
          return 0;
        else if ((elem == 1)|(elem == 3))
          return pow(0.5,d);
      }
    case 8:
      if(g == 1){ // AA*
        if ((elem == 2)|(elem == 4))
          return pow(0.5,d);
        else if ((elem == 1)|(elem == 3))
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 1)|(elem == 3))
          return 0;
        else if ((elem == 2)|(elem == 4))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if ((elem == 1)|(elem == 3))
          return 0;
        else if ((elem == 2)|(elem == 4))
          return pow(0.5,d);
      }
    } // end of Switch
  }
  return -1;
}


// Function for extracting entries of the emission probability matrix
// when the OPGP are considered the baseline (and so phase is unknown and the r.f's are sex-specific)
double Qentry_up(int config,int g,int d,int elem){
  // Rprintf("OPGP = %i; g = %i; d = %i; elem = %i\n",OPGP,g,d,elem);
  // Check first whether the genotype is missing. If so, return 1
  if (g == 4)
    return 1;
  else{ // If not missing proceed to compute value
    switch(config){
    case 1: // Informative
      if(g == 1){ // AA*
        if(elem == 1) 
          return 0;
        else if ((elem == 2)|(elem == 3))
          return pow(0.5,d);
        else if (elem == 4)
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 1)|(elem == 4))
          return 0;
        else if ((elem == 2)|(elem == 3))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if (elem == 4)
          return 0;
        else if ((elem == 2)|(elem == 3))
          return pow(0.5,d);
        else if (elem == 1)
          return 1 ;
      }
    case 2: // Paternal segregating
      if(g == 1){ // AA*
        if ((elem == 1)|(elem == 2))
          return pow(0.5,d);
        else if ((elem == 3)|(elem == 4))
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 3)|(elem == 4))
          return 0;
        else if ((elem == 1)|(elem == 2))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if ((elem == 3)|(elem == 4))
          return 0;
        else if ((elem == 1)|(elem == 2))
          return pow(0.5,d);
      }
    case 3: // Maternal segregating
      if(g == 1){ // AA*
        if ((elem == 1)|(elem == 3))
          return pow(0.5,d);
        else if ((elem == 2)|(elem == 4))
          return 1;
      }
      else if (g == 2){ // AB*
        if ((elem == 2)|(elem == 4))
          return 0;
        else if ((elem == 1)|(elem == 3))
          return 1-pow(0.5,d-1);
      }
      else if (g == 3){ // BB*
        if ((elem == 2)|(elem == 4))
          return 0;
        else if ((elem == 1)|(elem == 3))
          return pow(0.5,d);
      }
    }
  }
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

//// likelihood 1:
// Not sex-specific (assumed equal) 
// r.f constrainted to range [0,1/2].
// OPGP's (or phase) are assumed to be known
SEXP ll_fs_c(SEXP r, SEXP epsilon, SEXP delta, SEXP genon, SEXP depth, SEXP OPGP, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int row, col, ind, snp, nInd_c, nSnps_c;
  double *pll, *pr, *pgenon, *pdepth, *pOPGP, delta_c, epsilon_c;
  double alpha[4], alphaTemp[4], tempMat[4][4], csum, tsum;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  delta_c = REAL(delta)[0];
  epsilon_c = REAL(epsilon)[0];
  // Define the pointers to the other input R variables
  pOPGP = REAL(OPGP);
  pgenon = REAL(genon);
  pdepth = REAL(depth);
  pr = REAL(r);  
  // Define the output variable
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  pll = REAL(ll);
  double llval = 0;
  
  // Now compute the likelihood
  for(ind = 0; ind < nInd_c; ind++){
    // Compute alpha0
    csum = 0;
    for(col = 0; col < 4; col++){
      alpha[col] = Qentry(pOPGP[0], pgenon[ind], pdepth[ind], col+1, epsilon_c, delta_c);
      csum = csum + alpha[col];
    }
    llval = llval + log(csum);
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps_c; snp++){
      // multiple the transition matrix and the emission matrix together 
      for(row = 0; row < 4; row++){
        for(col = 0; col < 4; col++){
          tempMat[row][col] = Tmat(row, col, pr[snp-1]) * Qentry(pOPGP[snp], pgenon[ind + nInd_c*snp], pdepth[ind + nInd_c*snp], col+1, epsilon_c, delta_c);
        }
      }
      // Multiply alpha by the temp matrix
      csum = 0;
      for(col = 0; col < 4; col++){
        tsum = 0;
        for(row = 0; row < 4; row++)
          tsum = tsum + alpha[row]*tempMat[row][col];
        alphaTemp[col] = tsum;
        // Calculate the contribution to likelihood and add to cumulative lkielihood value
        csum = csum + alphaTemp[col];
      }
      llval = llval + log(csum);
      // copy the element of alphaTemp to alpha
      for(row = 0; row < 4; row++){
        alpha[row] = alphaTemp[row];
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
SEXP ll_fs_ss_c(SEXP r, SEXP genon, SEXP depth, SEXP OPGP, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int row, col, ind, snp, nInd_c, nSnps_c;
  double *pll, *pr, *pgenon, *pdepth, *pOPGP;
  double alpha[4], alphaTemp[4], tempMat[4][4], csum, tsum;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = REAL(OPGP);
  pgenon = REAL(genon);
  pdepth = REAL(depth);
  pr = REAL(r);  
  // Define the output variable
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  pll = REAL(ll);
  double llval = 0;
  
  // Now compute the likelihood
  for(ind = 0; ind < nInd_c; ind++){
    // Compute alpha0
    csum = 0;
    for(col = 0; col < 4; col++){
      alpha[col] = Qentry(pOPGP[0], pgenon[ind], pdepth[ind], col+1);
      csum = csum + alpha[col];
    }
    llval = llval + log(csum);
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps_c; snp++){
      // multiple the transition matrix and the emission matrix together 
      for(row = 0; row < 4; row++){
        for(col = 0; col < 4; col++){
          tempMat[row][col] = Tmat_ss(row, col, pr[snp-1], pr[snp-1+nSnps_c-1]) * Qentry(pOPGP[snp], pgenon[ind + nInd_c*snp], pdepth[ind + nInd_c*snp], col+1);
        }
      }
      // Multiply alpha by the temp matrix
      csum = 0;
      for(col = 0; col < 4; col++){
        tsum = 0;
        for(row = 0; row < 4; row++)
          tsum = tsum + alpha[row]*tempMat[row][col];
        alphaTemp[col] = tsum;
        // Calculate the contribution to likelihood and add to cumulative lkielihood value
        csum = csum + alphaTemp[col];
      }
      llval = llval + log(csum);
      // copy the element of alphaTemp to alpha
      for(row = 0; row < 4; row++){
        alpha[row] = alphaTemp[row];
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
SEXP ll_fs_up_ss_c(SEXP r, SEXP genon, SEXP depth, SEXP config, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int row, col, ind, snp, nInd_c, nSnps_c;
  double *pll, *pr, *pgenon, *pdepth, *pconfig;
  double alpha[4], alphaTemp[4], tempMat[4][4], csum, tsum;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pconfig = REAL(config);
  pgenon = REAL(genon);
  pdepth = REAL(depth);
  pr = REAL(r);  
  // Define the output variable
  SEXP ll;
  PROTECT(ll = allocVector(REALSXP, 1));
  pll = REAL(ll);
  double llval = 0;
  
  // Now compute the likelihood
  for(ind = 0; ind < nInd_c; ind++){
    // Compute alpha0
    csum = 0;
    for(col = 0; col < 4; col++){
      alpha[col] = Qentry_up(pconfig[0], pgenon[ind], pdepth[ind], col+1);
      csum = csum + alpha[col];
    }
    llval = llval + log(csum);
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps_c; snp++){
      // multiple the transition matrix and the emission matrix together 
      for(row = 0; row < 4; row++){
        for(col = 0; col < 4; col++){
          tempMat[row][col] = Tmat_ss(row, col, pr[snp-1], pr[snp-1+nSnps_c-1]) * Qentry_up(pconfig[snp], pgenon[ind + nInd_c*snp], pdepth[ind + nInd_c*snp], col+1);
        }
      }
      // Multiply alpha by the temp matrix
      csum = 0;
      for(col = 0; col < 4; col++){
        tsum = 0;
        for(row = 0; row < 4; row++)
          tsum = tsum + alpha[row]*tempMat[row][col];
        alphaTemp[col] = tsum;
        // Calculate the contribution to likelihood and add to cumulative lkielihood value
        csum = csum + alphaTemp[col];
      }
      llval = llval + log(csum);
      // copy the element of alphaTemp to alpha
      for(row = 0; row < 4; row++){
        alpha[row] = alphaTemp[row];
      }
    }
  }
  
  pll[0] = -1*llval;
  // Clean up and return likelihood value
  UNPROTECT(1);
  return ll;
}


//// Scaled versions of the likelihood to deal with overflow issues


//// likelihood 1:
// Not sex-specific (assumed equal) 
// r.f constrainted to range [0,1/2].
// OPGP's (or phase) are assumed to be known
SEXP ll_fs_scaled_c(SEXP r, SEXP genon, SEXP depth, SEXP OPGP, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int s1, s2, ind, snp, nInd_c, nSnps_c;
  double *pll, *pr, *pgenon, *pdepth, *pOPGP;
  double alphaTilde[4], alphaDot[4], sum, w_logcumsum, w_new;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = REAL(OPGP);
  pgenon = REAL(genon);
  pdepth = REAL(depth);
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
      alphaDot[s1] = 0.25 * Qentry(pOPGP[0], pgenon[ind], pdepth[ind], s1+1);
      sum = sum + alphaDot[s1];
    }
    // Scale forward probabilities
    for(s1 = 0; s1 < 4; s1++){
      alphaTilde[s1] = alphaDot[s1]/sum;
    }

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
        alphaDot[s2] = Qentry(pOPGP[snp], pgenon[ind + nInd_c*snp], pdepth[ind + nInd_c*snp], s2+1) * sum;
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
SEXP ll_fs_ss_scaled_c(SEXP r, SEXP genon, SEXP depth, SEXP OPGP, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int s1, s2, ind, snp, nInd_c, nSnps_c;
  double *pll, *pr, *pgenon, *pdepth, *pOPGP;
  double alphaTilde[4], alphaDot[4], sum, w_logcumsum, w_new;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = REAL(OPGP);
  pgenon = REAL(genon);
  pdepth = REAL(depth);
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
      alphaDot[s1] = 0.25 * Qentry(pOPGP[0], pgenon[ind], pdepth[ind], s1+1);
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
        alphaDot[s2] = Qentry(pOPGP[snp], pgenon[ind + nInd_c*snp], pdepth[ind + nInd_c*snp], s2+1) * sum;
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
SEXP ll_fs_up_ss_scaled_c(SEXP r, SEXP genon, SEXP depth, SEXP OPGP, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int s1, s2, ind, snp, nInd_c, nSnps_c;
  double *pll, *pr, *pgenon, *pdepth, *pOPGP;
  double alphaTilde[4], alphaDot[4], sum, w_logcumsum, w_new;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = REAL(OPGP);
  pgenon = REAL(genon);
  pdepth = REAL(depth);
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
      alphaDot[s1] = 0.25 * Qentry_up(pOPGP[0], pgenon[ind], pdepth[ind], s1+1);
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
        alphaDot[s2] = Qentry_up(pOPGP[snp], pgenon[ind + nInd_c*snp], pdepth[ind + nInd_c*snp], s2+1) * sum;
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


