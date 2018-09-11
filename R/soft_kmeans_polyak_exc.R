softkmeanspol_src <- '
using namespace arma;
using namespace Rcpp;

// Load all of the function inputs
mat data_a = as<mat>(data_r);
mat mean_hold = as<mat>(mean_r);
int maxit_a = as<int>(maxit_r);
double rate_a = as<double>(rate_r);
double base_a = as<double>(base_r);
int batch_a = as<int>(batch_r);
int groups_a = as<int>(groups_r);

// Get necessary dimensional elements
int obs_a = data_a.n_cols;
int dim_a = data_a.n_rows;

// Memory control
rowvec pi_a = ones<rowvec>(groups_a);
pi_a = pi_a/groups_a;
mat mean_a = mean_hold;
cube cov_a = zeros<cube>(dim_a,dim_a,groups_a);
for (int gg = 0; gg < groups_a; gg++) {
cov_a.slice(gg) = eye<mat>(dim_a,dim_a);
}

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

// Initialize gain
double gain = 1;

// Create polyak variables
rowvec pol_pi = pi_a;
mat pol_mean = mean_a;
cube pol_cov = cov_a;

// Initialize tau matrix
mat tau = zeros<mat>(groups_a,batch_a);

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
  
  // Convert back to regular parameters
  for (int gg = 0; gg < groups_a; gg++) {
    mean_a.col(gg) = T2.col(gg)/T1(gg);
  }
  
  // Compute polyak averages
  pol_mean = pol_mean*count/(count+1) + mean_a/(count+1);

  // Reset the model parameters
  model.set_means(mean_a);

}

// Initialize the Gaussian mixture model object with polyak components
gmm_full model_pol;
model_pol.reset(dim_a,groups_a);

// Set to polyak model parameters
model_pol.set_hefts(pol_pi);
model_pol.set_means(pol_mean);
model_pol.set_fcovs(pol_cov);

// Outputs
return Rcpp::List::create(
  Rcpp::Named("reg_means")=model.means,
  Rcpp::Named("pol_means")=model_pol.means);
'

softkmeans_pol <- cxxfunction(signature(data_r='numeric',
                                  mean_r='numeric',
                                  maxit_r='integer',
                                  groups_r='integer',
                                  rate_r='numeric',
                                  base_r='numeric',
                                  batch_r='integer'),
                              softkmeanspol_src, plugin = 'RcppArmadillo')

softkmean_clust_src <- '
using namespace arma;
using namespace Rcpp;

// Load all of the function inputs
mat data_a = as<mat>(data_r);
mat mean_a = as<mat>(mean_r);

// Get necessary dimensional elements
int dim_a = data_a.n_rows;
int groups_a = mean_a.n_cols;

// Make pi and cov
rowvec pi_a = ones<rowvec>(groups_a);
pi_a = pi_a/groups_a;
cube cov_a = zeros<cube>(dim_a,dim_a,groups_a);
for (int gg = 0; gg < groups_a; gg++) {
cov_a.slice(gg) = eye<mat>(dim_a,dim_a);
}

// Initialize the Gaussian mixture model object
gmm_full model;
model.reset(dim_a,groups_a);

// Set the model parameters
model.set_params(mean_a, cov_a, pi_a);

// Outputs
return Rcpp::List::create(
Rcpp::Named("allocations")=model.assign(data_a, eucl_dist));
'

softkmeans_clust <- cxxfunction(signature(data_r='numeric',
                                        mean_r='numeric'),
                              softkmean_clust_src, plugin = 'RcppArmadillo')

softkmean_ss_src <- '
using namespace arma;
using namespace Rcpp;

// Load all of the function inputs
mat data_a = as<mat>(data_r);
mat mean_a = as<mat>(mean_r);
uvec alloc_a = as<uvec>(alloc_r);

// Get necessary dimensional elements
int obs_a = data_a.n_cols;
int groups_a = mean_a.n_cols;

// Initialize an SS vector;
vec ss_a = zeros<vec>(groups_a);

// Do the computation loop
for (int nn = 0; nn < obs_a; nn++) {
for (int gg = 0; gg < groups_a; gg++) {
if (alloc_a(nn) == gg) {
ss_a(gg) = ss_a(gg) + 1 + pow(norm(data_a.col(nn)-mean_a.col(gg)),2);
}
}
}

// Outputs
return Rcpp::List::create(
Rcpp::Named("Within_SS")=ss_a);
'

softkmeans_ss <- cxxfunction(signature(data_r='numeric',
                                          mean_r='numeric',
                                       alloc_r='integer'),
                                softkmean_ss_src, plugin = 'RcppArmadillo')

Soft <- softkmeans_pol(t(Data), msEst$parameters$mean,
         1000,5,0.6,1,1000)

Clust <- softkmeans_clust(t(Data),Soft$reg_means)

SS_calc <- softkmeans_ss(t(Data),Soft$reg_means,Clust$allocations)
