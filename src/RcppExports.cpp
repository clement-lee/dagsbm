// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// reorder_dense
const arma::mat reorder_dense(const arma::mat Y, const arma::uvec sigma);
RcppExport SEXP _dagsbm_reorder_dense(SEXP YSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(reorder_dense(Y, sigma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dagsbm_reorder_dense", (DL_FUNC) &_dagsbm_reorder_dense, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_dagsbm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}