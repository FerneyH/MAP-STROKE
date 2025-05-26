#include <RcppArmadillo.h>
#include "cmath"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
//// NOT CURRENTLY USED
arma::mat armaExpCpp(arma::mat x){
  return(exp(x));
}

// [[Rcpp::export]]
//// NOT CURRENTLY USED
arma::mat armaExpCpp2(NumericMatrix x){
  arma::mat X(x.begin(), x.rows(), x.cols(), false, true);
  return(exp(X));
}

// [[Rcpp::export]]
//// CURRENTLY USED
//// Directly modifies a matrix x with the exponential function applied to each element
double armaExpCpp3(NumericMatrix x){
  arma::mat X(x.begin(), x.rows(), x.cols(), false, true);
  X = exp(X);
  return(0);
}


// [[Rcpp::export]]
//// CURRENTLY USED
//// Multiplies each column in X by the corresponding element in bycol
arma::mat colwiseMultCpp(NumericMatrix x, arma::vec bycol){
  arma::mat X(x.begin(), x.rows(), x.cols());
  for(unsigned int i = 0; i < X.n_cols; i++){
    X.col(i) *= bycol(i);
  }
  return(X);
}
