/*
##########################################################################
# Genotyping Uncertainty with Sequencing data and linkage MAPping (GUSMap)
# Copyright 2017-2020 Timothy P. Bilton <timothy.bilton@agresearch.co.nz>
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

#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
//#include "probFun.h"


/////////////////////////////////////////////////////////////////////////
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

 
///////////////////////////////////////////////////////////////

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
  double tempV, thres = -1e20;
  //Rcpp::Rcout << " - thres :" << thres << std::endl;
  for(ind = 0; ind < nInd; ind++){
    // Compute forward probabilities at snp 1
    //sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      tempV = log(Qentry(OPGP[0], pAA(ind, 0), pAB(ind, 0), pBB(ind, 0), s1+1));
      if(tempV - thres )
        eta(s1, 0) = log(0.25) + tempV;
      else
        eta(s1, 0) = thres;
      //if(ind < 2){
      //  Rcpp::Rcout << "ind :" << ind  << " - eta :" << eta(s1, 0) << std::endl;
      //  Rcpp::Rcout << "ind :" << ind  << " - log2 :" << tempV << std::endl;
      //}
      //sum = sum + alphaDot[s1];
    }
    
    for(snp = 1; snp < nSnps; snp++){
      // compute the next forward probabilities for snp \ell
      for(s2 = 0; s2 < 4; s2++){
        for(s1 = 0; s1 < 4; s1++){
          etaTemp[s1] = eta(s1, snp - 1) + log(Tmat(s1, s2, rf[snp-1]));
          if(etaTemp[s1] < thres)
            etaTemp[s1] = thres;
          //if(ind < 2){
          //  Rcpp::Rcout << "snp :" << snp << " - Tmat :" << Tmat(s1, s2, rf[snp-1]) << std::endl;
          //  Rcpp::Rcout << "snp :" << snp << " - eta :" << eta(s1, snp) << std::endl;
          //}
        }
        eta(s2, snp) = etaTemp[which_max(etaTemp)] + log(Qentry(OPGP[snp], pAA(ind, snp), pAB(ind, snp), pBB(ind, snp), s2+1));
        if(eta(s2, snp) < thres)
          eta(s2, snp) = thres;
        //if(ind < 2){
        //  Rcpp::Rcout << "snp :" << snp << " - eta :" << eta(s2, snp) << std::endl;
        //}
      }
    }
    
    for(s1 = 0; s1 < 4; s1++){
      etaTemp[s1] = eta(s1, nSnps - 1);
      //Rcpp::Rcout << "snp :" << snp << " - log(eta) :" << log(eta(s1, nSnps - 1)) << std::endl;
      if(etaTemp[s1] < thres)
        etaTemp[s1] = thres;
    }
    
    inferState = which_max(etaTemp);
    //Rcpp::Rcout << "snp :" << snp << " - etaTemp :" << which_max(etaTemp) << std::endl;
    states(ind, nSnps - 1) = inferState;
    // Now compute the most likely states in reverse
    for(snp = nSnps - 2; snp > -1; snp--){
      for(s1 = 0; s1 < 4; s1++){
        etaTemp[s1] = eta(s1, snp) + log(Tmat(s1, inferState, rf[snp-1]));
        if(etaTemp[s1] < thres)
          etaTemp[s1] = thres;
      }
      inferState = which_max(etaTemp);
      states(ind, snp) = inferState;
    }
  }
  return states;
}
  
  
  
