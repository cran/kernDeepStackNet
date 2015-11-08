// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List crossprodRcpp(Eigen::Map<Eigen::MatrixXd> A) {
  const int n(A.cols());
  Eigen::MatrixXd AtA(Eigen::MatrixXd(n, n).setZero().
                        selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint()));
  return Rcpp::List::create(Rcpp::Named("crossprod(A)") = AtA);
}
