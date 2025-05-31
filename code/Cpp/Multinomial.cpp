#include <RcppArmadillo.h>
#include "cmath"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
//// CURRENTLY USED (called by multinomialCubeCpp)
//// This is performs the inverse multinomial logit transformation on a matrix x 
//// that has already been exponentiated
arma::mat multilogitinverseCpp(arma::mat X){
  arma::mat res(X.n_rows, X.n_cols + 1);
  arma::vec temp = 1 / (1 + sum(X, 1));
  res.cols(1, res.n_cols - 1) = X.each_col() % temp;
  res.col(0) = temp;
  return(res);
}

// [[Rcpp::export]]
//// CURRENTLY USED
//// This function performs the inverse multinomial logit transformation on a cube linpreds
arma::cube multinomialCubeCpp(arma::cube linpreds){
  linpreds = exp(linpreds);
  arma::cube res(linpreds.n_slices, linpreds.n_rows, linpreds.n_cols + 1);
#ifdef _OPENMP
#pragma omp parallel for num_threads(8)
#endif
  for(unsigned int i = 0; i < linpreds.n_slices; i++){
    res.row(i) = multilogitinverseCpp(linpreds.slice(i));
  }
  return(res);
} 

// [[Rcpp::export]]
//// NOT CURRENTLY USED
//// Performs the same operations as mulinomialcubeCpp, but tries to reduce copies
arma::cube multinomialCubeCpp2(arma::cube linpreds){
  linpreds = exp(linpreds);
  arma::cube res(linpreds.n_slices, linpreds.n_rows, linpreds.n_cols + 1);
#ifdef _OPENMP
#pragma omp parallel for num_threads(8)
#endif
  for(unsigned int i = 0; i < linpreds.n_slices; i++){
    arma::vec temp = 1 / (1 + sum(linpreds.slice(i), 1));
    res.subcube(i, 0, 1, i, res.n_cols - 1, res.n_slices - 1) = linpreds.slice(i).each_col() % temp;
    res.subcube(i, 0, 0, i, res.n_cols - 1, 0) = temp;
  }
  return(res);
}

// [[Rcpp::export]]
// CURRENTLY USED
//// This function performs tcrossprod(beta, x) separately for each slice of beta in parallel
arma::cube arrayMultCpp(arma::cube beta, arma::mat x){
  arma::mat xt = x.t();
  arma::cube results(x.n_rows, beta.n_slices, beta.n_rows);
  //// Can change dimension of beta to be (2, 1, 3) to avoid transposes, but doesn't help much
  // arma::cube results(x.n_rows, beta.n_slices, beta.n_cols);
  // Might need to modify nthreads if beta.n_slices ever increases past 6
  unsigned int nthreads = beta.n_slices;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads)
#endif
  for(unsigned int i = 0; i < beta.n_slices; i++){
    results.col(i) = (beta.slice(i) * xt).t();
    // results.col(i) = x * beta.slice(i);
  }
  return(results);
}

// [[Rcpp::export]]
//// CURRENTLY USED
//// This function calculate the mean of the given cube x along the 0th dimension
arma::mat cubeMeanCpp(arma::cube x){
  return(arma::mean(x, 0));
}
