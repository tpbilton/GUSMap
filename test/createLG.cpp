//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// // [[Rcpp::export]]
double mean_top_n(NumericVector y, int n, IntegerVector rows) {
  NumericVector x = clone(y);
  x = x[rows];
  // sort x in ascending order
  std::sort(x.begin(), x.end());
  NumericVector top_n = tail(x,n); 
  return mean(top_n);
}

NumericVector find_maxLOD(IntegerVector unmapped, IntegerVector LG, NumericMatrix LOD, int nComp){
  NumericVector out(2);
  NumericVector meanLOD(unmapped.size());
  int nSnps = unmapped.size();
  int col;
  for(col=0; col < nSnps; col++){
    meanLOD[col] = mean_top_n(LOD(_,unmapped(col)), nComp, LG);
    Rcpp::Rcout << "meanLOD :" << meanLOD[col] << std::endl;
  }
  out[0] = which_max(meanLOD);
  out[1] = meanLOD[out[0]];
  return out;
}

// [[Rcpp::export]]
Rcpp::LogicalVector logical_index(Rcpp::IntegerVector idx, R_xlen_t n) {
  bool invert = false; 
  Rcpp::LogicalVector result(n, false);
  
  for (R_xlen_t i = 0; i < idx.size(); i++) {
    if (!invert && idx[i] < 0) invert = true;
    result[std::abs(idx[i])] = true;
  }
  
  if (!invert) return !result;
  return !result;
}



// [[Rcpp::export]]
List createLG_cpp(IntegerVector unmapped, NumericMatrix LOD, double LODthres, int nComp){
  
  // create the list to return
  List LG = List::create();
  
  // run algorithm
  //bool finish = false;
  double maxLOD, newLOD, maxMeanLOD;
  int row, col, rowSNP, colSNP, snp, nSnps, newSNP;
  IntegerVector newLG, dropSNPs;
  // initialize the new LG by finding the pair with the largest LOD score
  rowSNP = unmapped(1);
  colSNP = unmapped(0);
  maxLOD = LOD(rowSNP,colSNP);
  for(row = 2; row < unmapped.size(); row++){
    for(col = 0; col < row; col++){
      newLOD = LOD(unmapped(row),unmapped(col));
      if(newLOD > maxLOD){
        rowSNP = row;
        colSNP = col;
        maxLOD = newLOD;
      }
    }
  }
  Rcpp::Rcout << "maxLOD :" << maxLOD << std::endl;
  while(maxLOD > LODthres){
    // form the new LG 
    newLG = IntegerVector::create(unmapped(colSNP),unmapped(rowSNP));
    Rcpp::Rcout << "newLG :" << newLG << std::endl;
    // Drop these SNPs from the unmapped vector
    dropSNPs = IntegerVector::create(rowSNP,colSNP);
    unmapped = unmapped[logical_index(dropSNPs, unmapped.size())];
    Rcpp::Rcout << "unmapped :" << unmapped << std::endl;
    
    // find the next highet LOD score
    nSnps = unmapped.size();
    if(nSnps > 0){
      NumericVector meanLOD(nSnps);
      for(snp=0; snp < nSnps; snp++){
        meanLOD[snp] = mean_top_n(LOD(_,unmapped(snp)), nComp, newLG);
      }
      newSNP = which_max(meanLOD);
      Rcpp::Rcout << "newSNP :" << newSNP << std::endl;
      maxMeanLOD = meanLOD[newSNP];
      Rcpp::Rcout << "maxMeanLOD :" << maxMeanLOD << std::endl;
      // see whether we add the SNP
      while(maxMeanLOD > LODthres){
        // Add the new SNP to the LG
        newLG.push_back(unmapped[newSNP]);
        Rcpp::Rcout << "newLG :" << newLG << std::endl;
        // Drop the new SNP from the list of unmapped SNPs
        dropSNPs = IntegerVector::create(newSNP);
        unmapped = unmapped[logical_index(dropSNPs, unmapped.size())];
        Rcpp::Rcout << "unmapped :" << unmapped << std::endl;
        // if there are still unmapped SNPs, find the new highest mean LOD value
        if(unmapped.size() > 0){
          nSnps = unmapped.size();
          NumericVector meanLOD(nSnps);
          for(snp=0; snp < nSnps; snp++){
            meanLOD[snp] = mean_top_n(LOD(_,unmapped(snp)), nComp, newLG);
          }
          newSNP = which_max(meanLOD);
          Rcpp::Rcout << "newSNP :" << newSNP << std::endl;
          maxMeanLOD = meanLOD[newSNP];
          Rcpp::Rcout << "maxMeanLOD :" << maxMeanLOD << std::endl;
        } else{
          maxMeanLOD = -999;
        }
      }
    }
    // Find the highest LOD between the unmapped SNPs
    if(unmapped.size() > 1){
      maxLOD = -1;
      for(row = 1; row < unmapped.size(); row++){
        for(col = 0; col < row; col++){
          newLOD = LOD(unmapped(row),unmapped(col));
          if(newLOD > maxLOD){
            rowSNP = row;
            colSNP = col;
            maxLOD = newLOD;
          }
        }
      }
    }
    else
      maxLOD = -999;
    LG.push_back(newLG);
  }
  return LG;
}

