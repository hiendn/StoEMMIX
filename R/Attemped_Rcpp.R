## Load in C libraries
library(Rcpp)
library(RcppArmadillo)
library(inline)

stoEMMIX_src <- '

// Read data
Rcpp::NumericMatrix Xr(Xs);
Rcpp::NumericVector Cov(Covs);

// Get the number of rows and columns
int n = Xr.nrow();
int d = Xr.ncol();

// Create arma objects
arma::cube Cov_a(Cov)

// Return results
return Rcpp::List::create(Rcpp::Named("number of rows") = n,
                          Rcpp::Named("number of columns") = d,
                          Rcpp::Named("Covariance matricies") = Cov_a);

'

Simple_Code <- cxxfunction(signature(Xs='numeric',
                                     Covs='numeric'),
                           stoEMMIX_src,
                           plugin = 'RcppArmadillo')