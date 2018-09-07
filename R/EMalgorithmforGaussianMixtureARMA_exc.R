gmm_full_src <- '
using namespace arma;

// Convert necessary matrix objects to arma
mat data_a = as<mat>(data_r);
rowvec pi_a = as<rowvec>(pi_r);
mat mean_a = as<mat>(mean_r);
cube cov_a = as<cube>(cov_r);
int maxit_a = as<int>(maxit_r);
int groups_a = as<int>(groups_r);

// Initialize a gmm_full object
gmm_full model;

// Set the parameters
model.set_params(mean_a, cov_a, pi_a);
model.learn(data_a, groups_a, maha_dist, keep_existing, 0, maxit_a, 2.2e-16, false);

//
return Rcpp::List::create(
Rcpp::Named("log-likelihood")=model.sum_log_p(data_a),
Rcpp::Named("proportions")=model.hefts,
Rcpp::Named("means")=model.means,
Rcpp::Named("covariances")=model.fcovs);
'

GMM_arma <- cxxfunction(signature(data_r='numeric',
                                  pi_r='numeric',
                                  mean_r='numeric',
                                  cov_r='numeric',
                                  maxit_r='integer',
                                  groups_r='integer'),
                        gmm_full_src, plugin = 'RcppArmadillo')

GMM_arma(t(Data), msEst$parameters$pro, msEst$parameters$mean,
         msEst$parameters$variance$sigma,
         20,5)