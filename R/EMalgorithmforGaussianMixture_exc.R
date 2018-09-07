EMMIX_src <- '
// Convert all of the matrix objects to C objects
Rcpp::NumericMatrix data_c(data_r);
Rcpp::NumericVector pi_c(pi_r);
Rcpp::NumericMatrix mean_c(mean_r);
Rcpp::NumericVector cov_c(cov_r);

// Extract and convert the integer objects
int maxit_c(maxit_r);
int groups_c(groups_r);
int n_c = data_c.nrow();
int d_c = data_c.ncol();

// Convert necessary matrix objects to arma
arma::mat data_a(data_c.begin(),n_c,d_c,false);
arma::colvec pi_a(pi_c.begin(),groups_c,false);
arma::mat mean_a(mean_c.begin(),d_c,groups_c,false);

// Start the loop
for (int count = 0; count < maxit_c, count++){

// Compute Tau
tau_a = arma::mat(n_rows, n_cols, fill::zeros)
for (int obs = 0; obs < n_c, obs++) {
for (int comp = 0; comp < groups_c, comp++) {

// Extract the current covariance that I need
arma::mat arma

tau_a(obs,comp) = pi_c

}
}
}
'