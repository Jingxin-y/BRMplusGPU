// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//' @importFrom Rcpp evalCpp
//' @useDynLib brm
//' @exportPattern ˆ[[:alpha:]]+
// [[Rcpp::export]]
arma::mat mat_vec_mul(const arma::mat & m, const arma::colvec & v) {
  arma::mat r = m * v;
  return r;
}