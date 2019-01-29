# Load libaries
library(Rcpp)
library(RcppArmadillo)
library(inline)

# Inline C Source for Minibatch EM algorithm without truncation
stoEMMIXpol_src <- '
using namespace arma;
using namespace Rcpp;

// Load all of the function inputs
mat data_a = as<mat>(data_r);
rowvec pi_hold = as<rowvec>(pi_r);
mat mean_hold = as<mat>(mean_r);
cube cov_hold = as<cube>(cov_r);
int maxit_a = as<int>(maxit_r);
double rate_a = as<double>(rate_r);
double base_a = as<double>(base_r);
int batch_a = as<int>(batch_r);
int groups_a = as<int>(groups_r);

// Memory control
rowvec pi_a = pi_hold;
mat mean_a = mean_hold;
cube cov_a = cov_hold;

// Get necessary dimensional elements
int obs_a = data_a.n_cols;
int dim_a = data_a.n_rows;

// Initialize the Gaussian mixture model object
gmm_full model;
model.reset(dim_a,groups_a);

// Set the model parameters
model.set_params(mean_a, cov_a, pi_a);

// Initialize the sufficient statistics
rowvec T1 = batch_a*pi_a;
mat T2 = zeros<mat>(dim_a,groups_a);
for (int gg = 0; gg < groups_a; gg++) {
  T2.col(gg) = mean_a.col(gg)*pi_a(gg);
}
cube T3 = zeros<cube>(dim_a,dim_a,groups_a);
for (int gg = 0; gg < groups_a; gg++) {
  T3.slice(gg) = T1(gg)*cov_a.slice(gg) + T2.col(gg)*trans(T2.col(gg))/T1(gg);
}

// Initialize gain
double gain = 1;

// Create polyak variables
rowvec pol_pi = pi_a;
mat pol_mean = mean_a;
cube pol_cov = cov_a;

// Initialize tau matrix
mat tau = zeros<mat>(groups_a,batch_a);

// Loading bar setup
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

// Begin loop
for (int count = 0; count < maxit_a; count++) {
  
  // Update gain function
  gain = base_a*pow(count+1,-rate_a);
  
  // Construct a sample from the data
  IntegerVector seq_c = sample(obs_a,batch_a,0);
  uvec seq_a = as<uvec>(seq_c) - 1;
  mat subdata_a = data_a.cols(seq_a);
  
  // Compute the tau scores for the subsample
  for (int gg = 0; gg < groups_a; gg++) {
    tau.row(gg) = pi_a(gg)*exp(model.log_p(subdata_a,gg));
  }
  for (int nn = 0; nn < batch_a; nn++) {
    tau.col(nn) = tau.col(nn)/sum(tau.col(nn));
  }
  
  // Compute the new value of T1
  T1 = (1-gain)*T1 + gain*trans(sum(tau,1));
  
  // Compute the new value of T2
  for (int gg = 0; gg < groups_a; gg++) {
    T2.col(gg) = (1-gain)*T2.col(gg);
    for (int nn = 0; nn < batch_a; nn++) {
      T2.col(gg) = T2.col(gg) + gain*tau(gg,nn)*subdata_a.col(nn);
    }
  }
  
  // Compute the new value of T3
  for (int gg = 0; gg < groups_a; gg++) {
    T3.slice(gg) = (1-gain)*T3.slice(gg);
    for (int nn = 0; nn < batch_a; nn++) {
      T3.slice(gg) = T3.slice(gg) + gain*tau(gg,nn)*subdata_a.col(nn)*trans(subdata_a.col(nn));
    }
  }
  
  // Convert back to regular parameters
  pi_a = T1/batch_a;
  for (int gg = 0; gg < groups_a; gg++) {
    mean_a.col(gg) = T2.col(gg)/T1(gg);
    cov_a.slice(gg) = (T3.slice(gg)-T2.col(gg)*trans(T2.col(gg))/T1(gg))/T1(gg);
  }
  
  // Compute polyak averages
  pol_pi = pol_pi*count/(count+1) + pi_a/(count+1);
  pol_mean = pol_mean*count/(count+1) + mean_a/(count+1);
  pol_cov = pol_cov*count/(count+1) + cov_a/(count+1);
  
  // Reset the model parameters
  model.set_hefts(pi_a);
  model.set_means(mean_a);
  model.set_fcovs(cov_a);
  
  // Loading bar
  double percentage = (count+1)/maxit_a;
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);

}

// Initialize the Gaussian mixture model object with polyak components
gmm_full model_pol;
model_pol.reset(dim_a,groups_a);

// Set to polyak model parameters
model_pol.set_hefts(pol_pi);
model_pol.set_means(pol_mean);
model_pol.set_fcovs(pol_cov);

return Rcpp::List::create(
  Rcpp::Named("reg_log-likelihood")=model.sum_log_p(data_a),
  Rcpp::Named("reg_proportions")=model.hefts,
  Rcpp::Named("reg_means")=model.means,
  Rcpp::Named("reg_covariances")=model.fcovs,
  Rcpp::Named("pol_log-likelihood")=model_pol.sum_log_p(data_a),
  Rcpp::Named("pol_proportions")=model_pol.hefts,
  Rcpp::Named("pol_means")=model_pol.means,
  Rcpp::Named("pol_covariances")=model_pol.fcovs);
'

# Inline C source for Minibatch EM algorithm with truncation
stoEMMIXpoltrunc_src <- '
using namespace arma;
using namespace Rcpp;

// Load all of the function inputs
mat data_a = as<mat>(data_r);
rowvec pi_hold = as<rowvec>(pi_r);
mat mean_hold = as<mat>(mean_r);
cube cov_hold = as<cube>(cov_r);
int maxit_a = as<int>(maxit_r);
double rate_a = as<double>(rate_r);
double base_a = as<double>(base_r);
int batch_a = as<int>(batch_r);
int groups_a = as<int>(groups_r);
double c1_a = as<double>(c1_r);
double c2_a = as<double>(c2_r);
double c3_a = as<double>(c3_r);

// Memory control
rowvec pi_a = pi_hold;
mat mean_a = mean_hold;
cube cov_a = cov_hold;

// Get necessary dimensional elements
int obs_a = data_a.n_cols;
int dim_a = data_a.n_rows;

// Initialize the Gaussian mixture model object
gmm_full model;
model.reset(dim_a,groups_a);

// Set the model parameters
model.set_params(mean_a, cov_a, pi_a);

// Initialize the sufficient statistics
rowvec T1 = batch_a*pi_a;
rowvec T1_old_a = batch_a*pi_a;
mat T2 = zeros<mat>(dim_a,groups_a);
mat T2_old_a = zeros<mat>(dim_a,groups_a);
for (int gg = 0; gg < groups_a; gg++) {
  T2.col(gg) = mean_a.col(gg)*pi_a(gg);
  T2_old_a.col(gg) = mean_a.col(gg)*pi_a(gg);
}
cube T3 = zeros<cube>(dim_a,dim_a,groups_a);
cube T3_old_a = zeros<cube>(dim_a,dim_a,groups_a);
for (int gg = 0; gg < groups_a; gg++) {
  T3.slice(gg) = T1(gg)*cov_a.slice(gg) + T2.col(gg)*trans(T2.col(gg))/T1(gg);
  T3_old_a.slice(gg) = T1(gg)*cov_a.slice(gg) + T2.col(gg)*trans(T2.col(gg))/T1(gg);
}

// Initialize gain
double gain = 1;

// Initialize m for truncation
double mm = 0;

// Create polyak variables
rowvec pol_pi = pi_a;
mat pol_mean = mean_a;
cube pol_cov = cov_a;

// Initialize tau matrix
mat tau = zeros<mat>(groups_a,batch_a);

// Loading bar setup
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

// Begin loop
for (int count = 0; count < maxit_a; count++) {
  
  // Update gain function
  gain = base_a*pow(count+1,-rate_a);
  
  // Construct a sample from the data
  IntegerVector seq_c = sample(obs_a,batch_a,0);
  uvec seq_a = as<uvec>(seq_c) - 1;
  mat subdata_a = data_a.cols(seq_a);
  
  // Compute the tau scores for the subsample
  for (int gg = 0; gg < groups_a; gg++) {
    tau.row(gg) = pi_a(gg)*exp(model.log_p(subdata_a,gg));
  }
  for (int nn = 0; nn < batch_a; nn++) {
    tau.col(nn) = tau.col(nn)/sum(tau.col(nn));
  }
  
  // Compute the new value of T1
  T1 = (1-gain)*T1 + gain*trans(sum(tau,1));
  
  // Compute the new value of T2
  for (int gg = 0; gg < groups_a; gg++) {
    T2.col(gg) = (1-gain)*T2.col(gg);
    for (int nn = 0; nn < batch_a; nn++) {
      T2.col(gg) = T2.col(gg) + gain*tau(gg,nn)*subdata_a.col(nn);
    }
  }
  
  // Compute the new value of T3
  for (int gg = 0; gg < groups_a; gg++) {
    T3.slice(gg) = (1-gain)*T3.slice(gg);
    for (int nn = 0; nn < batch_a; nn++) {
      T3.slice(gg) = T3.slice(gg) + gain*tau(gg,nn)*subdata_a.col(nn)*trans(subdata_a.col(nn));
    }
  }
  
  // Truncation
  // Compute the minimum value of pi_a
  double pi_min_a = pi_a.min();
  // Compute the minimum/max value of mean_a
  double mean_abs_a = abs(mean_a).max();
  // Compute the minimum and maximum eigenvalues of cov_a
  double eigen_max_a = c3_a + 1000;
  double eigen_min_a = 0;
  for (int gg = 0; gg < groups_a; gg++) {
    vec eigen_val_a = eig_sym(cov_a.slice(gg));
    double internal_min_a = eigen_val_a.min();
    double internal_max_a = eigen_val_a.max();
    if (internal_min_a<eigen_min_a) {
      eigen_min_a = internal_min_a;
    }
    if (internal_max_a>eigen_max_a) {
      eigen_max_a = internal_max_a;
    }
  }
  
  // Perform truncation if necessary
  int do_truncate_a = (pi_min_a<(1/(c1_a+mm)))*(mean_abs_a>(c2_a+mm))*(eigen_min_a<(1/(c3_a+mm)))*(eigen_min_a>(c3_a+mm));
  if (do_truncate_a==1) {
    T1 = T1_old_a;
    T2 = T2_old_a;
    T3 = T3_old_a;
    mm = mm + 1;
  }
  
  // Convert back to regular parameters
  pi_a = T1/batch_a;
  for (int gg = 0; gg < groups_a; gg++) {
    mean_a.col(gg) = T2.col(gg)/T1(gg);
    cov_a.slice(gg) = (T3.slice(gg)-T2.col(gg)*trans(T2.col(gg))/T1(gg))/T1(gg);
  }
  
  
  // Compute polyak averages
  pol_pi = pol_pi*count/(count+1) + pi_a/(count+1);
  pol_mean = pol_mean*count/(count+1) + mean_a/(count+1);
  pol_cov = pol_cov*count/(count+1) + cov_a/(count+1);
  
  // Reset the model parameters
  model.set_hefts(pi_a);
  model.set_means(mean_a);
  model.set_fcovs(cov_a);
  
  // Loading bar
  double percentage = (count+1)/maxit_a;
  int val = (int) (percentage * 100);
  int lpad = (int) (percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush (stdout);
}

// Initialize the Gaussian mixture model object with polyak components
gmm_full model_pol;
model_pol.reset(dim_a,groups_a);

// Set to polyak model parameters
model_pol.set_hefts(pol_pi);
model_pol.set_means(pol_mean);
model_pol.set_fcovs(pol_cov);

return Rcpp::List::create(
  Rcpp::Named("reg_log-likelihood")=model.sum_log_p(data_a),
  Rcpp::Named("reg_proportions")=model.hefts,
  Rcpp::Named("reg_means")=model.means,
  Rcpp::Named("reg_covariances")=model.fcovs,
  Rcpp::Named("pol_log-likelihood")=model_pol.sum_log_p(data_a),
  Rcpp::Named("pol_proportions")=model_pol.hefts,
  Rcpp::Named("pol_means")=model_pol.means,
  Rcpp::Named("pol_covariances")=model_pol.fcovs);
'

# Define R function for algorithm without truncation
stoEMMIX_pol <- cxxfunction(signature(data_r='numeric',
                                      pi_r='numeric',
                                      mean_r='numeric',
                                      cov_r='numeric',
                                      maxit_r='integer',
                                      groups_r='integer',
                                      rate_r='numeric',
                                      base_r='numeric',
                                      batch_r='integer'),
                            stoEMMIXpol_src, plugin = 'RcppArmadillo')

# Define R function for algorithm with truncation
stoEMMIX_poltrunc <- cxxfunction(signature(data_r='numeric',
                                           pi_r='numeric',
                                           mean_r='numeric',
                                           cov_r='numeric',
                                           maxit_r='integer',
                                           groups_r='integer',
                                           rate_r='numeric',
                                           base_r='numeric',
                                           batch_r='integer',
                                           c1_r='numeric',
                                           c2_r='numeric',
                                           c3_r='numeric'),
                              stoEMMIXpoltrunc_src, plugin = 'RcppArmadillo')