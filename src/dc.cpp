#define ARMA_DONT_USE_OPENMP
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

//' Reorder rows and columns of a dense matrix with the same ordering
//'
//' @param Y A dense square matrix
//' @param sigma An integer vector that is a permutation of 0, 1, ..., ncol(Y)-1
//' @export
// [[Rcpp::export]]
const arma::mat reorder_dense(const arma::mat Y, const arma::uvec sigma) {
  return Y.submat(sigma, sigma);
}
