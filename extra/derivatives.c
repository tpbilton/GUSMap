#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

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


//Functions for computing the derivatives of the recombination fractions
double derRF(double rf, int s1, int s2){
  if(s1 == s2)
    return 2*rf - 2;
  else{
    int sSum = s1 + s2*4;
    if ((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12))
      return 2*rf;
    else
      return 1 - 2 * rf;
  }
}

// Functions for computing the derivate for the error parameter
double derEpsilon_geno(double epsilon, double depth_Ref, double depth, int geno){
  if(geno == 2){
    if(depth_Ref == 0)
      return depth * pow(epsilon, depth-1);
    else if (depth_Ref == depth)
      return -depth_Ref * pow(1-epsilon, depth_Ref-1);
    else
      return pow(1-epsilon, depth_Ref-1) * pow(epsilon, depth - depth_Ref - 1) * ((depth - depth_Ref) - depth*epsilon);
  }
  else if(geno == 1)
    return 0;
  else if(geno == 0){
    if(depth_Ref == 0)
      return -depth * pow(1-epsilon, depth-1);
    else if (depth_Ref == depth)
      return depth_Ref * pow(epsilon, depth_Ref-1);
    else
      return pow(1-epsilon, depth_Ref-1) * pow(epsilon, depth - depth_Ref - 1) * ((depth - depth_Ref) - depth*epsilon);
  }
  else
    return -9999;
}

double derEpsilon(double epsilon, double depth_Ref, double depth_Alt, double bin_coef, int OPGP, int elem){
  
  if( (depth_Ref + depth_Alt) == 0)
    return 0;
  else{
    switch(OPGP){
    case 1:
      if(elem == 1)
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 0);
      else if ((elem == 2)|(elem == 3))  
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
      else if (elem == 4)
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 2);
    case 2:
      if(elem == 3)
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 0);
      else if ((elem == 1)|(elem == 4))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
      else if (elem == 2)
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 2);
    case 3:
      if(elem == 2) 
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 0);
      else if ((elem == 1)|(elem == 4))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
      else if (elem == 3)
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 2);
    case 4:
      if(elem == 4) 
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 0);
      else if ((elem == 2)|(elem == 3))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
      else if (elem == 1)
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 2);
    case 5:
      if ((elem == 1)|(elem == 2))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
      else if ((elem == 3)|(elem == 4))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 2);
    case 6:
      if ((elem == 1)|(elem == 2))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 2);
      else if ((elem == 3)|(elem == 4))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
    case 7:
      if ((elem == 1)|(elem == 2))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 0);
      else if ((elem == 3)|(elem == 4))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
    case 8:
      if ((elem == 1)|(elem == 2))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
      else if ((elem == 3)|(elem == 4))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 0);
    case 9:
      if ((elem == 1)|(elem == 3))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
      else if ((elem == 2)|(elem == 4))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 2);
    case 10:
      if ((elem == 1)|(elem == 3))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 2);
      else if ((elem == 2)|(elem == 4))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
    case 11:
      if ((elem == 1)|(elem == 3))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 0);
      else if ((elem == 2)|(elem == 4))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
    case 12:
      if ((elem == 1)|(elem == 3))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
      else if ((elem == 2)|(elem == 4))
        return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 0);
    case 13:
      return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 2);
    case 14:
      return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
    case 15:
      return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 1);
    case 16:
      return bin_coef * derEpsilon_geno(epsilon, depth_Ref, depth_Ref+depth_Alt, 0);
    }
  }
  return -9999;
}


//// likelihood 3:
// Not sex-specific (assumed equal) 
// r.f constrainted to range [0,1/2].
// OPGP's (or phase) are assumed to be known
// Include error parameters
SEXP gradient_fs_scaled_err_c(SEXP r, SEXP epsilon, SEXP depth_Ref, SEXP depth_Alt, SEXP bin_coef,
                              SEXP Kaa, SEXP Kab, SEXP Kbb, SEXP OPGP, SEXP nInd, SEXP nSnps){
  // Initialize variables
  int s1, s2, ind, snp, nInd_c, nSnps_c, indx;
  double *pll, *pgrad, *pr, *pKaa, *pKab, *pKbb, *pOPGP, *pdepth_Ref, *pdepth_Alt, *pbin_coef;
  double alphaTilde[4], alphaDot[4], sum, w_logcumsum, w_new, epsilon_c, delta_temp;
  // Load R input variables into C
  nInd_c = INTEGER(nInd)[0];
  nSnps_c = INTEGER(nSnps)[0];
  // Define the pointers to the other input R variables
  pOPGP = REAL(OPGP);
  pKaa = REAL(Kaa);
  pKab = REAL(Kab);
  pKbb = REAL(Kbb);
  pdepth_Ref = REAL(depth_Ref);
  pdepth_Alt = REAL(depth_Alt);
  pbin_coef = REAL(bin_coef);
  pr = REAL(r);
  epsilon_c = REAL(epsilon)[0];

  // set up objects for the gradient
  int para;
  double phi[4][nSnps_c], phiDot[4][nSnps_c], sum_phi;
  for(para = 0; para < nSnps_c; para++){
    for(s1 = 0; s1 < 4; s1++)
      phi[s1][para] = 0;
  }
  
  // Define the output variable
  SEXP ll; 
  SEXP grad;
  ll = PROTECT(allocVector(REALSXP, 1));
  grad = PROTECT(allocVector(REALSXP, nSnps_c));
  pll = REAL(ll);
  pgrad = REAL(grad);
  double llval = 0;
  
  // Now compute the likelihood
  for(ind = 0; ind < nInd_c; ind++){
    // Compute forward probabilities at snp 1
    sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      //Rprintf("Q value :%.6f at snp %i in ind %i\n", Qentry(pOPGP[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1), 0, ind);
      alphaDot[s1] = 0.25 * Qentry(pOPGP[0], pKaa[ind], pKab[ind], pKbb[ind], s1+1);
      sum = sum + alphaDot[s1];
      // Compute the gradient for epsilon
      Rprintf("Derivative at snp %i and ind %i: %lf\n", 0, ind, derEpsilon(epsilon_c, pdepth_Ref[ind], pdepth_Alt[ind], pbin_coef[ind], pOPGP[0], s1+1));
      phi[s1][nSnps_c] = 0.25 * derEpsilon(epsilon_c, pdepth_Ref[ind], pdepth_Alt[ind], pbin_coef[ind], pOPGP[0], s1+1);
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
      indx = ind + nInd_c*snp;
      // compute the next forward probabilities for snp \ell
      for(s2 = 0; s2 < 4; s2++){
        delta_temp = Qentry(pOPGP[snp], pKaa[indx], pKab[indx], pKbb[indx], s2+1);
        sum = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum = sum + Tmat(s1, s2, pr[snp-1]) * alphaTilde[s1];
        }
        //Rprintf("Q value :%.6f at snp %i in ind %i\n", Qentry(pOPGP[snp], pKaa[ind + nInd_c*snp], pKab[ind + nInd_c*snp], pKbb[ind + nInd_c*snp], s2+1), snp, ind);
        alphaDot[s2] = delta_temp * sum;
      }
      // Compute the weight for snp \ell
      w_new = 0;
      for(s2 = 0; s2 < 4; s2++){
        w_new = w_new + alphaDot[s2];
      }
      // Add contribution to the likelihood
      llval = llval + log(w_new) + w_logcumsum;
      Rprintf("lik value: %lf\n", llval);
      w_logcumsum = w_logcumsum + log(w_new);
      
      // Compute the gradients for the parameters
      for(s2 = 0; s2 < 4; s2++){
        delta_temp = Qentry(pOPGP[snp], pKaa[indx], pKab[indx], pKbb[indx], s2+1);
        for(para = 0; para < nSnps_c; para++){
          sum_phi = 0;
          for(s1 = 0; s1 < 4; s1++){
            // gradient for r_{snp-1}
            if(para == snp - 1){
              sum_phi = alphaTilde[s1] * delta_temp * derRF(pr[snp-1], s1 , s2);
            }
            // gradient for epsilon
            else if(para == nSnps_c - 1){
              //Rprintf("Derivative at snp %i and ind %i for para %i: %lf\n", snp, ind, para, derEpsilon(epsilon_c, pdepth_Ref[indx], pdepth_Alt[indx], pbin_coef[indx], pOPGP[snp], s1+1));
              sum_phi = sum_phi + ( phi[s1][para]  + alphaTilde[s1] * 
                derEpsilon(epsilon_c, pdepth_Ref[indx], pdepth_Alt[indx], pbin_coef[indx], pOPGP[snp], s1+1) ) 
                * delta_temp * Tmat(s1, s2, pr[snp-1]);
            }
            // gradient for r_{k}, k=snp:nSnps_c
            else if(para < snp){
              //Rprintf("Derivative at snp %i and ind %i for para %i: %lf\n", snp, ind, para, phi[s1][para] * delta_temp * Tmat(s1, s2, pr[snp-1]));
              sum_phi = sum_phi + phi[s1][para] * delta_temp * Tmat(s1, s2, pr[snp-1]);
            }
          }
          Rprintf("sum_phi at snp %i and ind %i for para %i: %lf\n", snp, ind, para, sum_phi);
          phiDot[s2][para] = sum_phi;
        }
      }

      // Scale the forward probability vector
      for(s2 = 0; s2 < 4; s2++){
        alphaTilde[s2] = alphaDot[s2]/w_new;
        for(para = 0; para < nSnps_c; para++)
          phi[s2][para] = phiDot[s2][para]/exp(w_new);
      }
    }
    // Add the gradients on to
  }
  for(para = 0; para < nSnps_c; para++){
    sum_phi = 0;
    for(s1 = 0; s1 < 4; s1++){
      sum_phi = sum_phi + phi[s1][para];
    }
    pgrad[para] = sum_phi;
  }
  
  pll[0] = -1*llval;
  // Clean up and return likelihood value
  UNPROTECT(2);
  return grad;
}
