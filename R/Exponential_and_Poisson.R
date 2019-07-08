# Load libaries
library(Rcpp)
library(RcppArmadillo)
library(inline)

stoExponential_src <- '

using namespace arma;
using namespace Rcpp;

// Load all of the function inputs
NumericVector data_hold(data_r);
NumericVector pi_hold(pi_r);
NumericVector lambda_hold(lambda_r);
int maxit_a = as<int>(maxit_r);
double rate_a = as<double>(rate_r);
double base_a = as<double>(base_r);
int batch_a = as<int>(batch_r);
int groups_a = as<int>(groups_r);

// Get necessary dimensional elements
int obs_a = data_hold.length();

// Memory control
vec data_a = zeros<vec>(obs_a);
vec pi_a = zeros<vec>(groups_a);
vec lambda_a = zeros<vec>(groups_a);
data_a = data_hold;
pi_a = pi_hold;
lambda_a = lambda_hold;

// Initialize the sufficient statistics
vec T1 = batch_a*pi_a;
vec T2 = zeros<vec>(groups_a);
for (int gg = 0; gg < groups_a; gg++) {
  T2(gg) = T1(gg)/lambda_a(gg);
}

// Initialize gain
double gain = 1;

// Create polyak variables
vec pol_pi = pi_a;
vec pol_lambda = lambda_a;

// Initialize tau matrix
mat tau = zeros<mat>(groups_a,batch_a);

// Begin loop
for (int count = 0; count < maxit_a; count++) {
  
  // Update gain function
  gain = base_a*pow(count+1,-rate_a);
  
  // Construct a sample from the data
  IntegerVector seq_c = sample(obs_a,batch_a,0);
  uvec seq_a = as<uvec>(seq_c) - 1;
  vec subdata_a = data_a(seq_a);
  
  // Compute the tau scores for the subsample
  for (int gg = 0; gg < groups_a; gg++) {
    for (int nn = 0; nn < batch_a; nn++) {
    tau(gg,nn) = log(lambda_a(gg)*exp(-lambda_a(gg)*subdata_a(nn)));
    }
  }
  for (int nn = 0; nn < batch_a; nn++) {
    vec tau_current = tau.col(nn);
    double max_tau = tau_current.max();
    for (int gg = 0; gg < groups_a; gg++) {
      tau(gg,nn) = exp(log(pi_a(gg)) + tau(gg,nn)-max_tau);
    }
    tau.col(nn) = tau.col(nn)/sum(tau.col(nn));
  }
  
  // Compute the new value of T1
  for (int gg = 0; gg < groups_a; gg++) {
    T1(gg) = (1-gain)*T1(gg) + gain*sum(tau.row(gg));
  }
  
  // Compute the new value of T2
  for (int gg = 0; gg < groups_a; gg++) {
    T2(gg) = (1-gain)*T2(gg);
    for (int nn = 0; nn < batch_a; nn++) {
      T2(gg) = T2(gg) + gain*tau(gg,nn)*subdata_a(nn);
    }
  }
  
  // Convert back to regular parameters
  for (int gg = 0; gg < groups_a; gg++) {
    pi_a(gg) = T1(gg)/batch_a;
    lambda_a(gg) = T1(gg)/T2(gg);
  }
  
  // Compute polyak averages
  pol_pi = pol_pi*count/(count+1) + pi_a/(count+1);
  pol_lambda = pol_lambda*count/(count+1) + lambda_a/(count+1);
  
}


// Initialize the likelihoods
double log_likelihood = 0;
double log_likelihood_pol = 0;

// Initialize tau matrix
mat tau_like = zeros<mat>(groups_a,obs_a);
// Compute taus under regular parameters
  for (int gg = 0; gg < groups_a; gg++) {
    for (int nn = 0; nn < obs_a; nn++) {
    tau_like(gg,nn) = log(pi_a(gg)*lambda_a(gg)*exp(-lambda_a(gg)*data_a(nn)));
    }
  }
mat exp_tau = exp(tau_like);
rowvec log_sum_exp = log(sum(exp_tau,0));
log_likelihood = sum(log_sum_exp);

// Initialize tau matrix
tau_like = zeros<mat>(groups_a,obs_a);
// Compute taus under polyak parameters
  for (int gg = 0; gg < groups_a; gg++) {
    for (int nn = 0; nn < obs_a; nn++) {
    tau_like(gg,nn) = log(pol_pi(gg)*pol_lambda(gg)*exp(-pol_lambda(gg)*data_a(nn)));
    }
  }
exp_tau = exp(tau_like);
log_sum_exp = log(sum(exp_tau,0));
log_likelihood_pol = sum(log_sum_exp);

return Rcpp::List::create(
  Rcpp::Named("reg_log-likelihood")=log_likelihood,
  Rcpp::Named("reg_proportions")=pi_a,
  Rcpp::Named("reg_lambda")=lambda_a,
  Rcpp::Named("pol_log-likelihood")=log_likelihood_pol,
  Rcpp::Named("pol_proportions")=pol_pi,
  Rcpp::Named("pol_lambda")=pol_lambda);
'

EMExponential_src <- '

using namespace arma;
using namespace Rcpp;

// Load all of the function inputs
NumericVector data_hold(data_r);
NumericVector pi_hold(pi_r);
NumericVector lambda_hold(lambda_r);
int maxit_a = as<int>(maxit_r);
int groups_a = as<int>(groups_r);

// Get necessary dimensional elements
int obs_a = data_hold.length();

// Memory control
vec data_a = zeros<vec>(obs_a);
vec pi_a = zeros<vec>(groups_a);
vec lambda_a = zeros<vec>(groups_a);
data_a = data_hold;
pi_a = pi_hold;
lambda_a = lambda_hold;

// Initialize the sufficient statistics
vec T1 = obs_a*pi_a;
vec T2 = zeros<vec>(groups_a);
for (int gg = 0; gg < groups_a; gg++) {
  T2(gg) = T1(gg)/lambda_a(gg);
}

// Initialize tau matrix
mat tau = zeros<mat>(groups_a,obs_a);

// Begin loop
for (int count = 0; count < maxit_a; count++) {
  
  
  // Compute the tau scores for the subsample
  for (int gg = 0; gg < groups_a; gg++) {
    for (int nn = 0; nn < obs_a; nn++) {
      tau(gg,nn) = log(lambda_a(gg)*exp(-lambda_a(gg)*data_a(nn)));
    }
  }
  for (int nn = 0; nn < obs_a; nn++) {
    vec tau_current = tau.col(nn);
    double max_tau = tau_current.max();
    for (int gg = 0; gg < groups_a; gg++) {
      tau(gg,nn) = exp(log(pi_a(gg)) + tau(gg,nn)-max_tau);
    }
    tau.col(nn) = tau.col(nn)/sum(tau.col(nn));
  }
  
  // Compute the new value of T1
  for (int gg = 0; gg < groups_a; gg++) {
    T1(gg) = sum(tau.row(gg));
  }
  
  // Compute the new value of T2
  for (int gg = 0; gg < groups_a; gg++) {
    T2(gg) = 0;
    for (int nn = 0; nn < obs_a; nn++) {
      T2(gg) = T2(gg) + tau(gg,nn)*data_a(nn);
    }
  }
  
  // Convert back to regular parameters
  for (int gg = 0; gg < groups_a; gg++) {
    pi_a(gg) = T1(gg)/obs_a;
    lambda_a(gg) = T1(gg)/T2(gg);
  }
  
}


// Initialize the likelihoods
double log_likelihood = 0;

// Initialize tau matrix
mat tau_like = zeros<mat>(groups_a,obs_a);
// Compute taus under regular parameters
for (int gg = 0; gg < groups_a; gg++) {
  for (int nn = 0; nn < obs_a; nn++) {
    tau_like(gg,nn) = log(pi_a(gg)*lambda_a(gg)*exp(-lambda_a(gg)*data_a(nn)));
  }
}
mat exp_tau = exp(tau_like);
rowvec log_sum_exp = log(sum(exp_tau,0));
log_likelihood = sum(log_sum_exp);

return Rcpp::List::create(
  Rcpp::Named("reg_log-likelihood")=log_likelihood,
  Rcpp::Named("reg_proportions")=pi_a,
  Rcpp::Named("reg_lambda")=lambda_a);
'

stoPoisson_src <- '

using namespace arma;
using namespace Rcpp;

// Load all of the function inputs
NumericVector data_hold(data_r);
NumericVector pi_hold(pi_r);
NumericVector lambda_hold(lambda_r);
int maxit_a = as<int>(maxit_r);
double rate_a = as<double>(rate_r);
double base_a = as<double>(base_r);
int batch_a = as<int>(batch_r);
int groups_a = as<int>(groups_r);

// Get necessary dimensional elements
int obs_a = data_hold.length();

// Memory control
vec data_a = zeros<vec>(obs_a);
vec pi_a = zeros<vec>(groups_a);
vec lambda_a = zeros<vec>(groups_a);
data_a = data_hold;
pi_a = pi_hold;
lambda_a = lambda_hold;

// Initialize the sufficient statistics
vec T1 = batch_a*pi_a;
vec T2 = zeros<vec>(groups_a);
for (int gg = 0; gg < groups_a; gg++) {
  T2(gg) = lambda_a(gg)*T1(gg);
}

// Initialize gain
double gain = 1;

// Create polyak variables
vec pol_pi = pi_a;
vec pol_lambda = lambda_a;

// Initialize tau matrix
mat tau = zeros<mat>(groups_a,batch_a);

// Begin loop
for (int count = 0; count < maxit_a; count++) {
  
  // Update gain function
  gain = base_a*pow(count+1,-rate_a);
  
  // Construct a sample from the data
  IntegerVector seq_c = sample(obs_a,batch_a,0);
  uvec seq_a = as<uvec>(seq_c) - 1;
  vec subdata_a = data_a(seq_a);
  
  // Compute the tau scores for the subsample
  for (int gg = 0; gg < groups_a; gg++) {
    for (int nn = 0; nn < batch_a; nn++) {
    tau(gg,nn) = log(pow(lambda_a(gg),subdata_a(nn))*exp(-lambda_a(gg)))-lgamma(subdata_a(nn)+1);
    }
  }
  for (int nn = 0; nn < batch_a; nn++) {
    vec tau_current = tau.col(nn);
    double max_tau = tau_current.max();
    for (int gg = 0; gg < groups_a; gg++) {
      tau(gg,nn) = exp(log(pi_a(gg)) + tau(gg,nn)-max_tau);
    }
    tau.col(nn) = tau.col(nn)/sum(tau.col(nn));
  }
  
  // Compute the new value of T1
  for (int gg = 0; gg < groups_a; gg++) {
    T1(gg) = (1-gain)*T1(gg) + gain*sum(tau.row(gg));
  }
  
  // Compute the new value of T2
  for (int gg = 0; gg < groups_a; gg++) {
    T2(gg) = (1-gain)*T2(gg);
    for (int nn = 0; nn < batch_a; nn++) {
      T2(gg) = T2(gg) + gain*tau(gg,nn)*subdata_a(nn);
    }
  }
  
  // Convert back to regular parameters
  for (int gg = 0; gg < groups_a; gg++) {
    pi_a(gg) = T1(gg)/batch_a;
    lambda_a(gg) = T2(gg)/T1(gg);
  }
  
  // Compute polyak averages
  pol_pi = pol_pi*count/(count+1) + pi_a/(count+1);
  pol_lambda = pol_lambda*count/(count+1) + lambda_a/(count+1);
  
}


// Initialize the likelihoods
double log_likelihood = 0;
double log_likelihood_pol = 0;

// Initialize tau matrix
mat tau_like = zeros<mat>(groups_a,obs_a);
// Compute taus under regular parameters
  for (int gg = 0; gg < groups_a; gg++) {
    for (int nn = 0; nn < obs_a; nn++) {
    tau_like(gg,nn) = log(pi_a(gg)*pow(lambda_a(gg),data_a(nn))*exp(-lambda_a(gg)))-lgamma(data_a(nn)+1);
    }
  }
mat exp_tau = exp(tau_like);
rowvec log_sum_exp = log(sum(exp_tau,0));
log_likelihood = sum(log_sum_exp);

// Initialize tau matrix
tau_like = zeros<mat>(groups_a,obs_a);
// Compute taus under polyak parameters
  for (int gg = 0; gg < groups_a; gg++) {
    for (int nn = 0; nn < obs_a; nn++) {
    tau_like(gg,nn) = log(pol_pi(gg)*pow(pol_lambda(gg),data_a(nn))*exp(-pol_lambda(gg)))-lgamma(data_a(nn)+1);
    }
  }
exp_tau = exp(tau_like);
log_sum_exp = log(sum(exp_tau,0));
log_likelihood_pol = sum(log_sum_exp);

return Rcpp::List::create(
  Rcpp::Named("reg_log-likelihood")=log_likelihood,
  Rcpp::Named("reg_proportions")=pi_a,
  Rcpp::Named("reg_lambda")=lambda_a,
  Rcpp::Named("pol_log-likelihood")=log_likelihood_pol,
  Rcpp::Named("pol_proportions")=pol_pi,
  Rcpp::Named("pol_lambda")=pol_lambda);
'

EMPoisson_src <- '

using namespace arma;
using namespace Rcpp;

// Load all of the function inputs
NumericVector data_hold(data_r);
NumericVector pi_hold(pi_r);
NumericVector lambda_hold(lambda_r);
int maxit_a = as<int>(maxit_r);
int groups_a = as<int>(groups_r);

// Get necessary dimensional elements
int obs_a = data_hold.length();

// Memory control
vec data_a = zeros<vec>(obs_a);
vec pi_a = zeros<vec>(groups_a);
vec lambda_a = zeros<vec>(groups_a);
data_a = data_hold;
pi_a = pi_hold;
lambda_a = lambda_hold;

// Initialize the sufficient statistics
vec T1 = obs_a*pi_a;
vec T2 = zeros<vec>(groups_a);
for (int gg = 0; gg < groups_a; gg++) {
  T2(gg) = T1(gg)*lambda_a(gg);
}

// Initialize tau matrix
mat tau = zeros<mat>(groups_a,obs_a);

// Begin loop
for (int count = 0; count < maxit_a; count++) {
  
  
  // Compute the tau scores for the subsample
  for (int gg = 0; gg < groups_a; gg++) {
    for (int nn = 0; nn < obs_a; nn++) {
      tau(gg,nn) = log(pow(lambda_a(gg),data_a(nn))*exp(-lambda_a(gg)))-lgamma(data_a(nn)+1);
    }
  }
  for (int nn = 0; nn < obs_a; nn++) {
    vec tau_current = tau.col(nn);
    double max_tau = tau_current.max();
    for (int gg = 0; gg < groups_a; gg++) {
      tau(gg,nn) = exp(log(pi_a(gg)) + tau(gg,nn)-max_tau);
    }
    tau.col(nn) = tau.col(nn)/sum(tau.col(nn));
  }
  
  // Compute the new value of T1
  for (int gg = 0; gg < groups_a; gg++) {
    T1(gg) = sum(tau.row(gg));
  }
  
  // Compute the new value of T2
  for (int gg = 0; gg < groups_a; gg++) {
    T2(gg) = 0;
    for (int nn = 0; nn < obs_a; nn++) {
      T2(gg) = T2(gg) + tau(gg,nn)*data_a(nn);
    }
  }
  
  // Convert back to regular parameters
  for (int gg = 0; gg < groups_a; gg++) {
    pi_a(gg) = T1(gg)/obs_a;
    lambda_a(gg) = T2(gg)/T1(gg);
  }
  
}


// Initialize the likelihoods
double log_likelihood = 0;

// Initialize tau matrix
mat tau_like = zeros<mat>(groups_a,obs_a);
// Compute taus under regular parameters
for (int gg = 0; gg < groups_a; gg++) {
  for (int nn = 0; nn < obs_a; nn++) {
    tau_like(gg,nn) = log(pi_a(gg)*pow(lambda_a(gg),data_a(nn))*exp(-lambda_a(gg)))-lgamma(data_a(nn)+1);
  }
}
mat exp_tau = exp(tau_like);
rowvec log_sum_exp = log(sum(exp_tau,0));
log_likelihood = sum(log_sum_exp);

return Rcpp::List::create(
  Rcpp::Named("reg_log-likelihood")=log_likelihood,
  Rcpp::Named("reg_proportions")=pi_a,
  Rcpp::Named("reg_lambda")=lambda_a);
'

# Define R function for algorithm without truncation
stoExponential_pol <- cxxfunction(signature(data_r='numeric',
                                      pi_r='numeric',
                                      lambda_r='numeric',
                                      maxit_r='integer',
                                      groups_r='integer',
                                      rate_r='numeric',
                                      base_r='numeric',
                                      batch_r='integer'),
                                  stoExponential_src, plugin = 'RcppArmadillo')

# Define R function for algorithm without truncation
EMExponential_pol <- cxxfunction(signature(data_r='numeric',
                                            pi_r='numeric',
                                            lambda_r='numeric',
                                            maxit_r='integer',
                                            groups_r='integer'),
                                  EMExponential_src, plugin = 'RcppArmadillo')


# Define R function for algorithm without truncation
stoPoisson_pol <- cxxfunction(signature(data_r='numeric',
                                            pi_r='numeric',
                                            lambda_r='numeric',
                                            maxit_r='integer',
                                            groups_r='integer',
                                            rate_r='numeric',
                                            base_r='numeric',
                                            batch_r='integer'),
                                  stoPoisson_src, plugin = 'RcppArmadillo')

# Define R function for algorithm without truncation
EMPoisson_pol <- cxxfunction(signature(data_r='numeric',
                                           pi_r='numeric',
                                           lambda_r='numeric',
                                           maxit_r='integer',
                                           groups_r='integer'),
                                 EMPoisson_src, plugin = 'RcppArmadillo')

# Exponential clustering  
Exp_clust <- function(data, pi_vec, lambda_vec) {
  tau <- matrix(NA,length(data),length(pi_vec))
  for (gg in 1:length(pi_vec)) {
    tau[,gg] <- log(pi_vec[gg]) + dexp(data,lambda_vec[gg],log=TRUE)
  }
  return(apply(tau,1,which.max))
}

# Poisson clustering  
Pois_clust <- function(data, pi_vec, lambda_vec) {
  tau <- matrix(NA,length(data),length(pi_vec))
  for (gg in 1:length(pi_vec)) {
    tau[,gg] <- log(pi_vec[gg]) + dpois(data,lambda_vec[gg],log=TRUE)
  }
  return(apply(tau,1,which.max))
}

SAMPLE <- c(rexp(20000),rexp(5000,3),rexp(5000,9))
Sto_mix <- stoExponential_pol(SAMPLE,pi_r=c(4/6,1/6,1/6),lambda_r=c(1,3,9),
                   maxit_r = 3000,groups_r = 3, 0.6,1-10^-10, batch_r = 1000)
EM_mix <- EMExponential_pol(SAMPLE,pi_r=c(4/6,1/6,1/6),lambda_r=c(1,3,9),
                            maxit_r = 1000,groups_r = 3)
print(Sto_mix)
print(EM_mix)
CLUSTER <- Exp_clust(SAMPLE,EM_mix$reg_proportions,EM_mix$reg_lambda)
table(CLUSTER)


SAMPLE <- c(rpois(20000,1),rpois(5000,3),rpois(5000,9))
Sto_mix <- stoPoisson_pol(SAMPLE,pi_r=c(4/6,1/6,1/6),lambda_r=c(1,3,9),
                              maxit_r = 3000,groups_r = 3, 0.6,1-10^-10, batch_r = 1000)
EM_mix <- EMPoisson_pol(SAMPLE,pi_r=c(4/6,1/6,1/6),lambda_r=c(1,3,9),
                            maxit_r = 1000,groups_r = 3)
print(Sto_mix)
print(EM_mix)
CLUSTER <- Pois_clust(SAMPLE,EM_mix$reg_proportions,EM_mix$reg_lambda)
table(CLUSTER)
