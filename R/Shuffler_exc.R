library(Rcpp)
library(RcppArmadillo)
library(inline)

shuffle_src <- '
using namespace arma;
using namespace Rcpp;

int num_a = as<int>(num_r);
IntegerVector seq_c = seq_len(num_a);
IntegerVector seq_a = sample(num_a,num_a,0);

return Rcpp::List::create(
  Rcpp::Named("sequence") = seq_a);
'
Shuffler <- cxxfunction(signature(num_r='numeric'),
                        shuffle_src,
                        plugin = 'RcppArmadillo')

matrix_shuffle_src <- '
using namespace arma;
using namespace Rcpp;

int p_a = as<int>(p_r);
int up_a = as<int>(up_r);
mat mat_a = as<mat>(mat_r);
IntegerVector seq_c = seq_len(p_a);
IntegerVector seq_a = sample(seq_c,up_a,0);
uvec seq_arma = as<uvec>(seq_a) - 1;

mat submat = mat_a.cols(seq_arma);

return Rcpp::List::create(
Rcpp::Named("sequence") = seq_a,
Rcpp::Named("submat") = submat);
'
Matrix_Shuffler <- cxxfunction(signature(p_r='numeric',
                                  up_r='numeric',
                                  mat_r='numeric'),
                        matrix_shuffle_src,
                        plugin = 'RcppArmadillo')