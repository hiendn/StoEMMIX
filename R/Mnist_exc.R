library(Rcpp)
library(RcppArmadillo)
library(inline)

gmm_full_cluster_src <- '
using namespace arma;

// Convert necessary matrix objects to arma
mat data_a = as<mat>(data_r);
rowvec pi_a = as<rowvec>(pi_r);
mat mean_a = as<mat>(mean_r);
cube cov_a = as<cube>(cov_r);

// Initialize a gmm_full object
gmm_full model;

// Set the parameters
model.set_params(mean_a, cov_a, pi_a);

// Return
return Rcpp::List::create(
Rcpp::Named("Cluster")=model.assign(data_a, prob_dist));
'

GMM_arma_cluster <- cxxfunction(signature(data_r='numeric',
                                          pi_r='numeric',
                                          mean_r='numeric',
                                          cov_r='numeric'),
                                gmm_full_cluster_src, plugin = 'RcppArmadillo')

stoEMMIXpolsafe_src <- '
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
double safe_a = as<double>(safe_r);
int batch_a = as<int>(batch_r);
int groups_a = as<int>(groups_r);

// Memory control
rowvec pi_a = pi_hold;
mat mean_a = mean_hold;
cube cov_a = cov_hold;

// Get necessary dimensional elements
int obs_a = data_a.n_cols;
int dim_a = data_a.n_rows;

// Safe variances
for (int gg = 0; gg < groups_a; gg++) {
if (prod(cov_a.slice(gg).diag())==0) {
cov_a.slice(gg) = safe_a*eye<mat>(dim_a,dim_a);
}
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
rowvec tau_max = zeros<rowvec>(batch_a);

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
tau.row(gg) = log(pi_a(gg)) + model.log_p(subdata_a,gg);
}
tau_max = max(tau,0);
for (int nn = 0; nn < batch_a; nn++) {
tau.col(nn) = tau.col(nn) - tau_max(nn) - log(sum(exp(tau.col(nn)-tau_max(nn))));
}
tau = exp(tau);

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
if (prod(cov_a.slice(gg).diag())==0) {
cov_a.slice(gg) = safe_a*eye<mat>(dim_a,dim_a);
}
}

// Compute polyak averages
pol_pi = pol_pi*count/(count+1) + pi_a/(count+1);
pol_mean = pol_mean*count/(count+1) + mean_a/(count+1);
pol_cov = pol_cov*count/(count+1) + cov_a/(count+1);

// Reset the model parameters
model.set_hefts(pi_a);
model.set_means(mean_a);
model.set_fcovs(cov_a);

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

stoEMMIX_polsafe <- cxxfunction(signature(data_r='numeric',
                                  pi_r='numeric',
                                  mean_r='numeric',
                                  cov_r='numeric',
                                  maxit_r='integer',
                                  groups_r='integer',
                                  rate_r='numeric',
                                  base_r='numeric',
                                  batch_r='integer',
                                  safe_r='numeric'),
                        stoEMMIXpolsafe_src, plugin = 'RcppArmadillo')




Data_hold <- train[,-785]
PCA <- prcomp(Data_hold)
Data_hold <- predict(PCA,Data_hold)
Data <- Data_hold[,1:200]

Samp <- sample(1:10,dim(Data)[1],replace = T)
id <- Samp
KM <- kmeans(Data,centers = 10,iter.max = 100, nstart = 10)
id <- KM$cluster
K <- max(id)
# estimate mixture parameters
Pi <- prop.table(tabulate(id))
Mu <- t(sapply(1:K, function(k){ colMeans(Data[id == k,]) }))
S <- sapply(1:K, function(k){ var(Data[id == k,]) })
dim(S) <- c(dim(Data)[2], dim(Data)[2], K)
Sto <- stoEMMIX_polsafe(t(Data), Pi, t(Mu),
                    S,
                    1000,K,0.6,1,10000,10)
Cluster <- GMM_arma_cluster(t(Data),Sto$reg_proportions,Sto$reg_means,Sto$reg_covariances)
adjustedRandIndex(Cluster$Cluster,train$y)
adjustedRandIndex(KM$cluster,train$y)