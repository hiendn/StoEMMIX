gmm_full_src <- '
// Convert all of the matrix objects to C objects
Rcpp::NumericMatrix data_c(data_r);
Rcpp::NumericVector pi_c(pi_r);
Rcpp::NumericMatrix mean_c(mean_r);
Rcpp::NumericVector cov_c(cov_r);

// Extract and convert the integer objects
int n_c = data_c.nrow();
int d_c = data_c.ncol();

// Convert necessary matrix objects to arma
arma::mat data_a(data_c.begin(), n_c, d_c, false);
arma::rowvec pi_a(pi_c.begin(), groups_r, false);
arma::mat mean_a(mean_c.begin(), d_c, groups_r, false);
arma::cube cov_a(cov_c.begin(), d_c, d_c, groups_r, false);

// Initialize a gmm_full object
arma::gmm_full model;

// Set the parameters
model.set_params(mean_a, cov_a, pi_a);
bool model.learn(data_a, groups_r, maha_dist, keep_existing, 0, maxit_r, 2.2e-16, false);
'

GMM_arma <- cxxfunction(signature(data_r='numeric',
                                  pi_r='numeric',
                                  mean_r='numeric',
                                  cov_r='numeric',
                                  maxit_r='integer',
                                  groups_r='integer'),
                        gmm_full_src, plugin = 'RcppArmadillo')