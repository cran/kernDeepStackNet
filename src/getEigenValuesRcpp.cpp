// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
/**
 * @title Using Eigen for eigenvalues
 * @author Dirk Eddelbuettel
 * @license GPL (>= 2)
 * @tags eigen matrix
 * @summary This example shows how to compute eigenvalues using Eigen
 *

 * [A previous post](../armadillo-eigenvalues) showed how to compute 
 * eigenvalues using the [Armadillo](http://arma.sf.net) library via RcppArmadillo.
 *
 * Here, we do the same using [Eigen](http://eigen.tuxfamily.org) and
 * the RcppEigen package.
 *
 */

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd getEigenValuesRcpp(Eigen::Map<Eigen::MatrixXd> M) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
    return es.eigenvalues();
}
