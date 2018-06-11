
#include <math.h>
#include <stdio.h>
#include <limits.h>
#include "probFun.h"


// Function for computing binomial coefficients
// taken  from this website "https://rosettacode.org/wiki/Evaluate_binomial_coefficients#C"
static unsigned long gcd_ui(unsigned long x, unsigned long y) {
  unsigned long t;
  if (y < x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  // y1 <- x0 % y0 ; x1 <- y0
  }
  return x;
}

unsigned long binomial(unsigned long a, unsigned long b) {
  unsigned long n, d, g, r = 1;
  n = a + b;
  if (a == 0) return 1;
  if (a == 1) return n;
  if (a >= n) return (a == n);
  if (a > n/2) a = n-a;
  for (d = 1; d <= a; d++) {
    if (r >= ULONG_MAX/n) {  // Possible overflow 
      unsigned long nr, dr;  // reduced numerator / denominator 
      g = gcd_ui(n, d);  nr = n/g;  dr = d/g;
      g = gcd_ui(r, dr);  r = r/g;  dr = dr/g;
      if (r >= ULONG_MAX/nr) return 0;  // Unavoidable overflow
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
