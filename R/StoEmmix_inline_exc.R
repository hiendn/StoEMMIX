library(Rcpp)
library(RcppArmadillo)
library(inline)

stoEMMIX_src <- '
using namespace arma;
using namespace Rcpp;

// Load all of the function inputs
mat data_a = as<mat>(data_r);
rowvec pi_a = as<rowvec>(pi_r);
mat mean_a = as<mat>(mean_r);
cube cov_a = as<cube>(cov_r);
int maxit_a = as<int>(maxit_r);
double rate_a = as<double>(rate_r);
int batch_a = as<int>(batch_r);
int groups_a = as<int>(groups_r);

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

// Initialize tau matrix
mat tau = zeros<mat>(groups_a,batch_a);

// Begin loop
for (int count = 0; count < maxit_a; count++) {

// Update gain function
gain = pow(count+1,-rate_a);

// Construct a sample from the data
IntegerVector seq_c = seq_len(obs_a);
IntegerVector seq_a = sample(seq_c,batch_a,0);
uvec seq_arma = as<uvec>(seq_a) - 1;
mat subdata_a = data_a.cols(seq_arma);

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

// Reset the model parameters
model.set_hefts(pi_a);
model.set_means(mean_a);
model.set_fcovs(cov_a);

}

return Rcpp::List::create(
Rcpp::Named("log-likelihood")=model.sum_log_p(data_a),
Rcpp::Named("proportions")=model.hefts,
Rcpp::Named("means")=model.means,
Rcpp::Named("covariances")=cov_a);
'

stoEMMIX <- cxxfunction(signature(data_r='numeric',
                                  pi_r='numeric',
                                  mean_r='numeric',
                                  cov_r='numeric',
                                  maxit_r='integer',
                                  groups_r='integer',
                                  rate_r='numeric',
                                  batch_r='integer'),
                        stoEMMIX_src, plugin = 'RcppArmadillo')

# stoEMMIX(t(Data), msEst$parameters$pro, msEst$parameters$mean,
#          msEst$parameters$variance$sigma,
#          10000,5,0.6,100)
