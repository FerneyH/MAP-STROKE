#include <RcppArmadillo.h>
#include "cmath"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// This section contains functions used to calculate mean for truncated Poisson
// [[Rcpp::export]]
//// NOT CURRENTLY USED
//// This function is the simplest implementation
arma::mat truncPoisMeanCpp(NumericMatrix x, unsigned int max){
  arma::mat X(x.begin(), x.rows(), x.cols());
  for(unsigned int i = 0; i < X.n_cols; i++){
    for(unsigned int j = 0; j < X.n_rows; j++){
      // Getting lambda
      double lambda = X(j, i);
      
      // Creating objects to store expectation and probabilities
      //// Getting probability of x = 0
      double curprob = exp(-lambda);
      double totprob = curprob;
      double expectation = 0;
      
      // Getting rest of the probabilities
      for(unsigned int k = 1; k <= max; k++){
        curprob *= (lambda / k);
        totprob += curprob;
        expectation += curprob * k;
      }
      
      // Rescaling expectation
      expectation /= totprob;
      X(j, i) = expectation;
    }
  }
  return(X);
}

// [[Rcpp::export]]
//// NOT CURRENTLY USED
//// This function is the same as the above implementation, but does it in parallel by column
arma::mat truncPoisMeanCpp2(NumericMatrix x, unsigned int max){
  arma::mat X(x.begin(), x.rows(), x.cols());
#ifdef _OPENMP
#pragma omp parallel for num_threads(8)
#endif
  for(unsigned int i = 0; i < X.n_cols; i++){
    for(unsigned int j = 0; j < X.n_rows; j++){
      // Getting lambda
      double lambda = X(j, i);
      
      // Creating objects to store expectation and probabilities
      //// Getting probability of x = 0
      double curprob = exp(-lambda);
      double totprob = curprob;
      double expectation = 0;
      
      // Getting rest of the probabilities
      for(unsigned int k = 1; k <= max; k++){
        curprob *= (lambda / k);
        totprob += curprob;
        expectation += curprob * k;
      }
      
      // Rescaling expectation
      expectation /= totprob;
      X(j, i) = expectation;
    }
  }
  return(X);
}

// [[Rcpp::export]]
//// CURRENTLY USED
//// This function does it in parallel, but also directly modifies the supplied x 
//// matrix with the means
double truncPoisMeanCpp3(NumericMatrix x, unsigned int max){
  arma::mat X(x.begin(), x.rows(), x.cols(), false, true);
#ifdef _OPENMP
#pragma omp parallel for num_threads(8)
#endif
  for(unsigned int i = 0; i < X.n_cols; i++){
    for(unsigned int j = 0; j < X.n_rows; j++){
      // Getting lambda
      double lambda = X(j, i);
      
      // Creating objects to store expectation and probabilities
      //// Getting probability of x = 0
      double curprob = exp(-lambda);
      double totprob = curprob;
      double expectation = 0;
      
      // Getting rest of the probabilities
      for(unsigned int k = 1; k <= max; k++){
        curprob *= (lambda / k);
        totprob += curprob;
        expectation += curprob * k;
      }
      
      // Rescaling expectation
      expectation /= totprob;
      X(j, i) = expectation;
    }
  }
  return(0);
}

// [[Rcpp::export]]
//// CURRENTLY USED
//// This function generates random numbers from a truncated poisson 
//// based on the untruncated Poisson lambda parameter
arma::mat rtruncPoisCpp(NumericMatrix x, unsigned int max){
  arma::mat X(x.begin(), x.rows(), x.cols());
  arma::mat unifs = arma::randu(X.n_rows, X.n_cols);
  for(unsigned int i = 0; i < X.n_cols; i++){
    for(unsigned int j = 0; j < X.n_rows; j++){
      // Getting lambda
      double lambda = X(j, i);
      double curunif = unifs(j, i);
      arma::vec probs(max + 1, arma::fill::zeros);
      
      // Creating objects to store expectation and probabilities
      //// Getting probability of x = 0
      double curprob = exp(-lambda);
      double totprob = curprob;
      probs(0) = totprob;
      
      // Getting rest of the probabilities
      for(unsigned int k = 1; k <= max; k++){
        curprob *= (lambda / k);
        totprob += curprob;
        probs(k) = totprob;
      }
      
      // Rescaling expectation
      probs /= totprob;
      
      // Getting sampled number
      unsigned int k = 0;
      while(k <= max && probs(k) < curunif){
        k++;
      }
      
      X(j, i) = k;
    }
  }
  return(X);
}
