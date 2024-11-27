#define ARMA_DONT_USE_OPENMP
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;
using namespace arma;





// 00) Prelim & generic functions
void nan_to_minus_infinity(double & x) {
  // turn NaN x to -Inf
  if (isnan(x)) {
    x = -INFINITY;
  }
}

const double lnan(const double l) {
  return (l != l) ? -INFINITY : l; // or isnan(l)
}

template <class T>
const double lpg(T v) {
  // log of product of gamma of v (vector or matrix)
  return accu(lgamma(v));
}

const double lr1() {
  return log(runif(1)[0]);
}

const bool update(double & par_curr,
                  const double par_prop,
                  double & lpost_curr,
                  const double lpost_prop,
                  double & s,
                  const int i,
                  const int burn,
                  const double factor = 10.0,
                  const double adj = 0.0) {
  // M-H update
  const bool accept_reject = lr1() < lpost_prop - lpost_curr + adj;
  par_curr = accept_reject ? par_prop : par_curr;
  lpost_curr = accept_reject ? lpost_prop : lpost_curr;
  if (i < burn) {
    s = sqrt(s * s + (accept_reject ? 3.0 : (-1.0)) * s * s / factor / sqrt(i + 1.0));
  }
  return accept_reject;
}

const double nC2(int n) {
  return n * (n - 1.0) / 2.0;
}





// 01) Rcpp f()s
const int sample_1(const IntegerVector seq) {
  // sample 1 value from seq
  return Rcpp::RcppArmadillo::sample(seq, 1, true, NumericVector::create())[0];
}

const int sample_w(const IntegerVector seq, const NumericVector weights) {
  // sample 1 value from seq with weights
  return Rcpp::RcppArmadillo::sample(seq, 1, true, weights)[0];
}

const IntegerVector sl1(const int n) {
  return seq_len(n) - 1;
}

const NumericVector tv(const double x) {
  return NumericVector::create(x);
}

const IntegerVector ti(const int x) {
  return IntegerVector::create(x);
}

const LogicalVector tl(const bool x) {
  return LogicalVector::create(x);
}





// 02) arma f()s - manipulating matrices & vectors
const vec Z_to_N(const uvec Z, const int K) {
  // from memberships to group sizes
  // built-in hist_c() neck to neck
  vec N(K);
  uvec ui;
  for (int i = 0; i < K; i++) {
    ui = find(Z == i);
    N[i] = ui.size();
  }
  return N;
}

const mat Z_to_E(const sp_mat Y, const uvec Z, const int K) {
  // from memberships to edges between groups
  mat E(K, K);
  sp_mat Y1;
  uvec ui, uj;
  for (int i = 0; i < K; i++) {
    ui = find(Z == i);
    for (int j = 0; j < K; j++) {
      uj = find(Z == j);
      Y1 = Y.cols(uj).t();
      Y1 = Y1.cols(ui); // what if ui or uj empty
      E(i, j) = accu(nonzeros(Y1));
    }
  }
  return E;
}

List edges(const rowvec v, const colvec w, const uvec Z, const int K, const int n) {
  // lengths of v, w & Z should all equal n
  rowvec r(K);
  colvec c(K);
  uvec u;
  for (int j = 0; j < K; j++) {
    u = find(Z == j);
    if (u.size() == 0) {
      r[j] = 0.0;
      c[j] = 0.0;
    } else {
      r[j] = sum(v(u));
      c[j] = sum(w(u));
    }
  }
  List L = List::create(Named("r") = r,
                        Named("c") = c);
  return L;
}

template <class T>
T reorder (T Y, const uvec sigma) {
  // reorder rows & columns of Y simultaneously
  return Y.submat(sigma, sigma);
}

//' Reorder rows and columns of a dense matrix with the same ordering
//'
//' @param Y A dense square matrix
//' @param sigma An integer vector that is a permutation of 0, 1, ..., ncol(Y)-1
//' @export
// [[Rcpp::export]]
const arma::mat reorder_dense(const arma::mat Y, const arma::uvec sigma) {
  return Y.submat(sigma, sigma);
}

//' Reorder rows and columns of a sparse matrix with the same ordering
//'
//' @param Y A sparse square matrix
//' @param sigma An integer vector that is a permutation of 0, 1, ..., ncol(Y)-1
//' @export
// [[Rcpp::export]]
const arma::sp_mat reorder_sparse(const arma::sp_mat Y, const arma::uvec sigma) {
  const sp_mat Yt = Y.t(), Y0 = Yt.cols(sigma).t();
  return Y0.cols(sigma);
}

const bool check_equal_seq(const uvec U, const int l) {
  // check if U is identical to (0, 1, ..., l-1)
  const uvec seq1 = regspace<uvec>(0, l-1);
  return approx_equal(seq1, U, "absdiff", 0);
  // slightly faster than two alternatives below:
  // convert IntegerVector to uvec & approx_equal()
  // wrap() two IntegerVectors & comparing
}

const mat Z_to_M(const uvec Z_star, const vec xi_star, const int K, const int n, const bool dag = true) {
  // cross-ref: Z_to_E()
  // Z_star already topo. sorted
  mat Z_mat = xi_star * xi_star.t(), M(K, K);
  if (dag) {
    Z_mat = trimatu(Z_mat, 1); // omit main diag.
  }
  else {
    Z_mat.diag().zeros();
  }
  uvec ui, uj;
  for (int i = 0; i < K; i++) {
    ui = find(Z_star == i);
    for (int j = 0; j < K; j++) {
      uj = find(Z_star == j);
      M(i, j) = accu(Z_mat.submat(ui, uj));
    }
  }
  return M;
}

void check_NEM(const sp_mat Y, const int K, const uvec Z, const uvec sigma, const vec xi, const int n, const vec N_curr, const mat E_curr, const mat M_curr, const bool dag = true) {
  vec N_check = Z_to_N(Z, K);
  mat E_check = Z_to_E(Y, Z, K),
    M_check = Z_to_M(Z(sigma), xi(sigma), K, n, dag);
  Rcout << "Max. N diff = " << max(abs(N_check - N_curr)) << endl;
  Rcout << "Max. E diff = " << max(abs(vectorise(E_check - E_curr))) << endl;
  Rcout << "Max. M diff = " << max(abs(vectorise(M_check - M_curr))) << endl;
}

const int check_sigma_phi(const uvec sigma, const uvec phi) {
  const bool check = approx_equal(sort_index(sigma), phi, "both", 0.0, 0.0);
  return (int) !check; // 0 means equal
}

const double lV(const double x, const int n) {
  // log(x * (x + 1) * ... * (x + n - 1))
  if (n < 1) {
    stop("lV: n has to be greater than or equal to 1.");
  }
  return sum(log(regspace(x, x + n - 1)));
}

const double cal_entropy(const vec N, const int K, const int n) {
  // entropy of a vector of sizes
  if (N.size() != K) {
    stop("cal_entropy: length of N is not equal to K.");
  }
  if (sum(N) != n) {
    stop("cal_entropy: sum of N is not equal to n.");
  }
  const vec props = N / (n + 0.0);
  return -dot(props, log(props)); //// 2022-05-02: have to make it +ve by adding -ve sign in front
}

void check_for_lpost(const sp_mat Y,
                     const uvec Z,
                     const int K,
                     const uvec sigma,
                     const bool dag = true) {
  // generic checks for collapsed sampler, finite or infinite regime
  const int n = Y.n_rows;
  if (Y.n_cols != n) {
    stop("check_for_lpost: Y has to be a square matrix.");
  }
  if (Z.size() != n || sigma.size() != n) {
    stop("check_for_lpost: lengths of Z & sigma have to equal # rows/columns of Y.");
  }
  if (K > n || K <= 0) {
    stop("check_for_lpost: K has to be between 1 & the number of nodes (inclusive).");
  }
  // b) Y, Z, & sigma separately
  if (any(nonzeros(Y) < 0.0)) { // requires dense
    stop("check_for_lpost: all elements of Y have to be non-negative (integers).");
  }
  const uvec Z0 = unique(Z);
  if (!check_equal_seq(Z0, K)) {
    stop("check_for_lpost: Z has to have K unique values {0, 1, ..., K-1}.");
  }
  uvec sigma0 = sigma;
  std::sort(sigma0.begin(), sigma0.end());
  if (!check_equal_seq(sigma0, n)) {
    stop("check_for_lpost: sigma has to be a permutation of {0, 1, ..., n-1)}.");
  }
  // c) Y & sigma jointly
  if (dag) { // no need to check if not dag model
    const sp_mat Yt = Y.t();
    sp_mat Y0 = Yt.cols(sigma).t();
    Y0 = Y0.cols(sigma);
    if (!Y0.is_trimatu()) {
      stop("check_for_lpost: sigma has to be such that reordered Y is upper triangular.");
    }
  }
}

const double lcross(const mat A, const mat B, const int j) {
  const double l = 0.0
    + (lpg(A.row(j)) - dot(A.row(j), log(B.row(j))))
    + (lpg(A.col(j)) - dot(A.col(j), log(B.col(j))))
    - (lgamma(A(j, j)) - A(j, j) * log(B(j, j)));
  return l;
}

const double lincr(const mat P, const mat Q, const rowvec rp, const colvec cp, const rowvec rq, const colvec cq, const int j) {
  mat R = P, S = Q;
  R.row(j) += rp;
  R.col(j) += cp;
  S.row(j) += rq;
  S.col(j) += cq;
  return lcross(R, S, j) - lcross(P, Q, j);
}

const double lpg_2r2c(const mat A, const double c, const int i, const int j) {
  // sum of lgamma of rows & cols i & j of A + c
  uvec ij(2);
  ij[0] = i;
  ij[1] = j;
  const double
    plus = lpg(A.rows(ij) + c) + lpg(A.cols(ij) + c),
    minus = lpg(A.submat(ij, ij) + c); // double-counted
  return plus - minus;
}

const double spdl_2r2c(const mat A, const mat B, const double c, const mat D, const int i, const int j) {
  // sum of rows & cols i & j of (A - B) % log(c + D)
  // spdl stands for sumproduct of diff & log
  const double plus =
    dot(A.row(i) - B.row(i), log(D.row(i) + c)) +
    dot(A.row(j) - B.row(j), log(D.row(j) + c)) +
    dot(A.col(i) - B.col(i), log(D.col(i) + c)) +
    dot(A.col(j) - B.col(j), log(D.col(j) + c)),
    minus =
    (A(i, i) - B(i, i)) * log(D(i, i) + c) +
    (A(i, j) - B(i, j)) * log(D(i, j) + c) +
    (A(j, i) - B(j, i)) * log(D(j, i) + c) +
    (A(j, j) - B(j, j)) * log(D(j, j) + c);
  return plus - minus;
}

const double nC2_weighted(const uvec Z, const vec xi, const int i, const bool dag = true) {
  // similar to nC2() for Z's update
  const uvec u = find(Z == i);
  const vec xi_partial = xi(u);
  return (dag ? 0.5 : 1.0) * (pow(sum(xi_partial), 2.0) - dot(xi_partial, xi_partial));
}





// 03) collapsed sampler, both regimes
const double lpost_par(const vec N,
                       const int K,
                       const int n,
                       const vec xi,
                       const int k,
                       const double gamma,
                       const double theta,
                       const double alpha,
                       const double a0,
                       const double b0,
                       const bool finite,
                       const double p_fin,
                       const double a1,
                       const double b1,
                       const double aa0,
                       const double ba0,
                       const double ab0,
                       const double bb0,
                       const double axi,
                       const double bxi,
                       const double ak,
                       const double bk) {
  double lpost;
  if ((finite && (k < K || gamma <= 0.0)) ||
      (finite && p_fin == 0.0) ||
      (!finite && (alpha < 0.0 || alpha >= 1.0 || theta < -alpha)) ||
      (!finite && p_fin == 1.0) || a0 <= 0.0 || b0 <= 0.0) {
    lpost = -INFINITY;
  }
  else {
    vec p(K), q(K - 1), r(n - 1);
    const double
      alpha0 = finite ? (- gamma) : alpha,
      theta0 = finite ? (k * gamma) : theta;
    for (int i = 0; i < K; i++) {
      p[i] = (N[i] == 1.0) ? 0.0 : sum(log(-alpha0 + regspace(1.0, N[i] - 1.0)));
    }
    q = theta0 + alpha0 * regspace(1.0, K - 1.0);
    r = theta0 + regspace(1.0, n - 1.0);
    const NumericVector xi_rcpp = wrap(xi);
    lpost = sum(p) + sum(log(q)) - sum(log(r)) +
      dgamma(tv(a0), aa0, 1.0 / ba0, true)[0] +
      dgamma(tv(b0), ab0, 1.0 / bb0, true)[0] +
      sum(dgamma(xi_rcpp, axi, 1.0 / bxi, true));
    if (finite) {
      const IntegerVector k0 = ti(k);
      const NumericVector v0 = tv(gamma);
      lpost +=
        dnbinom(k0, ak, bk, true)[0] - log(1.0 - pow(bk, ak)) +
        dgamma(v0, a1, 1.0 / b1, true)[0] +
        log(p_fin);
    }
    else { // infinite regime
      const NumericVector v0 = tv(theta0 + alpha);
      lpost +=
        dgamma(v0, a1, 1.0 / b1, true)[0] +
        log(1.0 - p_fin);
    }
  }
  return lnan(lpost);
}

const double lpost_all(const sp_mat Y,
                       const uvec Z,
                       const int K,
                       const uvec sigma,
                       const vec xi,
                       const int k,
                       const double gamma,
                       const double theta,
                       const double alpha,
                       const double a0,
                       const double b0,
                       const bool finite,
                       const double p_fin,
                       const double a1,
                       const double b1,
                       const double aa0,
                       const double ba0,
                       const double ab0,
                       const double bb0,
                       const double axi,
                       const double bxi,
                       const double ak,
                       const double bk,
                       const bool dag = true) {
  // log-joint posterior
  // (k/theta, alpha) are PY pars (both regimes)
  const int n = Y.n_rows;
  check_for_lpost(Y, Z, K, sigma, dag);
  if (a1 <= 0.0 || b1 <= 0.0 ||
      aa0 <= 0.0 || ba0 <= 0.0 ||
      ab0 <= 0.0 || bb0 <= 0.0 ||
      axi <= 0.0 || bxi <= 0.0 ||
      ak <= 0.0 || bk <= 0.0) {
    stop("lpost: hyperparameters a1, b1, aa0, ba0, ab0, bb0, axi, bxi, ak & bk have to be positive.");
  }
  const vec N = Z_to_N(Z, K),
    log_xi = log(xi);
  mat log_W(n, n, fill::zeros);
  log_W.each_col() += log_xi;
  log_W.each_row() += log_xi.t();
  const mat M = Z_to_M(Z(sigma), xi(sigma), K, n, dag), E = Z_to_E(Y, Z, K);
  double lpost = lpost_par(N, K, n, xi, k, gamma, theta, alpha, a0, b0, finite, p_fin, a1, b1, aa0, ba0, ab0, bb0, axi, bxi, ak, bk)
    - (lgamma(a0) - a0 * log(b0)) * K * K
    + lpg(E + a0) - accu((E + a0) % log(M + b0))
    + dot(vectorise(log_W), vectorise(Y));
  // checks for k, gamma, theta & alpha in lpost_par()
  return lnan(lpost);
}

const int mod(const int x, const int y) {
  return x - x / y * y; // truncation towards 0, not flooring!
}

//' Shift the position-th element in sigma by distance
//'
//' @param sigma A permutation of 0, 1, ..., length(sigma)-1
//' @param position integer between 0 and length(sigma)-1 inclusive
//' @param distance integer between -length(sigma) and length(sigma) inclusive
//' @export
// [[Rcpp::export]]
const IntegerVector ulam(const arma::uvec sigma, const int position, const int distance) {
  // sigma is a permutation of {0, 1, ..., n-1}
  const int n = sigma.size();
  // 0 <= position <= n-1, -n < distance < n
  IntegerVector u = wrap(sigma);
  u.attr("dim") = R_NilValue;
  int value = u[position],
    i = mod(position + distance + n, n);
  // trun towards 0 instead of flooring, hence + n
  u.erase(position);
  u.insert(u.begin() + i, value);
  return u;
}

//' MCMC sampler of DAG-SBM, with possible model selection between finite and infinite regime
//'
//' @param Y A sparse square matrix
//' @param sigma An integer vector that is a permutation of 0, 1, ..., ncol(Y)-1
//' @param scalars Data frame of 1 row with the following columns: K, seed, iter, thin, burn, freq, node, scan, L, pg, p_fin, mean_k, a_gamma, b_gamma, a_theta, b_theta, a_alpha, b_alpha, a1, b1, aa0, ba0, ab0, bb0, axi, bxi, ak, bk, dag
//' @export
// [[Rcpp::export]]
List gvs(const arma::sp_mat Y,
         arma::uvec sigma,
         DataFrame scalars) {
  // GVS sampler for DAG SBM
  int K = scalars["K"];
  const int
    seed = scalars["seed"],
    iter = scalars["iter"],
    thin = scalars["thin"],
    burn = scalars["burn"],
    freq = scalars["freq"],
    node = scalars["node"],
    scan = scalars["scan"],
    L = scalars["L"];
  const double
    pg = scalars["pg"],
    p_fin = scalars["p_fin"],
    mean_k = scalars["mean_k"],
    a_gamma = scalars["a_gamma"],
    b_gamma = scalars["b_gamma"],
    a_theta = scalars["a_theta"],
    b_theta = scalars["b_theta"],
    a_alpha = scalars["a_alpha"],
    b_alpha = scalars["b_alpha"],
    a1 = scalars["a1"],
    b1 = scalars["b1"],
    aa0 = scalars["aa0"],
    ba0 = scalars["ba0"],
    ab0 = scalars["ab0"],
    bb0 = scalars["bb0"],
    axi = scalars["axi"],
    bxi = scalars["bxi"],
    ak = scalars["ak"],
    bk = scalars["bk"];
  const bool dag = scalars["dag"];
  Rcout << "Seed = " << seed << endl;
  Rcout << "Initial K = " << K << endl;
  Rcout << "After burn-in of length " << burn << ", " << endl;
  Rcout << "save every " << thin << " iteration(s), " << endl;
  Rcout << "to achieve a sample of " << iter << " iteration(s)" << endl;
  Rcout << "Printing log-posterior every " << freq << " iteration(s)" << endl;
  Rcout << "Prior: C_{ij} ~ Gamma(a0, b0)" << endl;
  Rcout << "Prior: Pr(finite) = " << p_fin << endl;
  Rcout << "Prior: a0 ~ Gamma(" << aa0 << ", " << ba0 << ")" << endl;
  Rcout << "Prior: b0 ~ Gamma(" << ab0 << ", " << bb0 << ")" << endl;
  Rcout << "Prior: xi ~ iid Gamma(" << axi << ", " << bxi << ")" << endl;
  if (dag) {
    Rcout << "Model is for directed ACYCLIC graphs" << endl;
  }
  else {
    Rcout << "Model is for directed graphs" << endl;
    sigma = regspace<uvec>(0, Y.n_rows - 1); // & stays fixed
  }
  Rcout << endl;
  // 01) checks & initial save
  // objects w/ dim n (x n)
  const int n = Y.n_rows;
  const sp_mat Yt = Y.t();
  const mat YY(Y + Yt);
  const vec sY = sum(YY, 1);
  sp_mat Y0(n, n);
  uvec seqn = regspace<uvec>(0, n - 1),
    Z_star = seqn - seqn / K * K,
    phi = sort_index(sigma),
    Z = Z_star(phi),
    Z0(n), Z1(n),
    node_set, // changing size
    indices; // changing size
  NumericVector xi_rcpp = rgamma(n, axi, 1.0 / bxi);
  vec xi = as<vec>(xi_rcpp), xi_star = xi(sigma);
  IntegerVector seqn_rcpp = wrap(seqn), pq(2), nodes;
  // for updating sigma non-incrementally
  if (L <= 0 || L >= n) {
    stop("gvs: L has to be between 0 & n exclusive.");
  }
  bool modulo;
  ivec seqL = regspace<ivec>(-L, L);
  IntegerVector seqL_rcpp = wrap(seqL),
    d(n); // will be proposed sigma
  seqL_rcpp.erase(seqL_rcpp.begin() + L);
  uvec sigma_prop(n);
  // scalars
  bool finite;
  if (p_fin < 0.0 || p_fin > 1.0) {
    stop("gvs: p_fin has to be between 0.0 & 1.0 inclusive.");
  }
  else if (p_fin == 0.0) {
    finite = false; // and stays
  }
  else if (p_fin == 1.0) {
    finite = true; // and stays
  }
  else {
    finite = false; // initial value
  }
  int k = K + 1;
  double gamma = 1.0, alpha = 0.5, theta = 0.5,
    a0 = 1.0, b0 = 1.0,
    lpost_prop, lpost_curr, lpost_orig;
  auto lpost0 = [n, p_fin, a1, b1, aa0, ba0, ab0, bb0, axi, bxi, ak, bk](const vec N, const int K, const vec xi, const int k, const double gamma, const double theta, const double alpha, const double a0, const double b0, const bool finite) {
                  return lpost_par(N, K, n, xi, k, gamma, theta, alpha, a0, b0, finite, p_fin, a1, b1, aa0, ba0, ab0, bb0, axi, bxi, ak, bk);
  };
  auto lpost = [p_fin, a1, b1, aa0, ba0, ab0, bb0, axi, bxi, ak, bk, dag](const sp_mat Y, const uvec Z, const int K, const uvec sigma, const vec xi, const int k, const double gamma, const double theta, const double alpha, const double a0, const double b0, const bool finite) {
                 return lpost_all(Y, Z, K, sigma, xi, k, gamma, theta, alpha, a0, b0, finite, p_fin, a1, b1, aa0, ba0, ab0, bb0, axi, bxi, ak, bk, dag);
  };
  lpost_curr = lpost(Y, Z, K, sigma, xi, k, gamma, theta, alpha, a0, b0, finite);
  Rcout << "Iter 0: Log-posterior = " << lpost_curr << endl;
  // 02) initialisation: for updating
  // omnipresent ones & big auxiliary matrices
  vec N = Z_to_N(Z, K), N0, N1;
  mat E = Z_to_E(Y, Z, K), E_prop(K, K), E0, E1,
    M = Z_to_M(Z(sigma), xi_star, K, n, dag), M_prop(K, K), M0, M1,
    F(n, K), F0, F1,
    G(K, n), G0, G1,
    H(n, K), U(n, K), U0, U1,
    I(K, n), V(K, n), V0, V1;
  List rc; uvec u, v;
  rowvec Yrp(n), Yro(n); colvec Ycp(n), Yco(n);
  for (int p = 0; p < n; p++) {
    Yrp = Yt.col(p).t();
    Ycp = Y.col(p);
    rc = edges(Yrp, Ycp, Z, K, n);
    rowvec r = rc["r"];
    colvec c = rc["c"];
    F.row(p) = r;
    G.col(p) = c;
  }
  if (dag) {
    for (int q = 0; q < n; q++) {
      for (int i = 0; i < K; i++) {
        if (q == 0) {
          v = find(Z_star.tail(n-1-q) == i) + q + 1;
          U(q, i) = xi_star[q] * accu(xi_star(v));
          V(i, q) = 0.0;
        }
        else if (q == n-1) {
          u = find(Z_star.head(q) == i);
          U(q, i) = 0.0;
          V(i, q) = xi_star[q] * accu(xi_star(u));
        }
        else {
          u = find(Z_star.head(q) == i);
          v = find(Z_star.tail(n-1-q) == i) + q + 1;
          U(q, i) = xi_star[q] * accu(xi_star(v));
          V(i, q) = xi_star[q] * accu(xi_star(u));
        }
      } 
    }
  }
  else {
    for (int q = 0; q < n; q++) {
      for (int i = 0; i < K; i++) {
        u = find(Z_star == i);
        U(q, i) = xi_star[q] * accu(xi_star(u));
        if (Z_star[q] == i) {
          U(q, i) -= pow(xi_star[q], 2.0);
        }
      }
    }
    V = U.t();
  }
  H = U.rows(phi);
  I = V.cols(phi);
  double Eij, Eji, Mij, Mji,
    a, b, ldiff, r, lA0, lA1,
    lpost_inf, lpost_fin, xip,
    s_gamma = 0.01, s_alpha = 0.01, s_theta = 0.01,
    s_a0 = 0.01, s_b0 = 0.01;
  vec s_xi(n), xi_prop(n);
  s_xi.fill(0.01);
  mat P(K, K), Q(K, K);
  rowvec rp(K), rq(K), rn(n);
  colvec cp(K), cq(K), cn(n);
  bool add_group, accept_reject;
  rowvec xir; colvec xic;
  // 03) initialisation: for saving
  umat Z_mat(iter, n), phi_mat(iter, n);
  mat xi_mat(iter, n);
  double lpost_max = lpost_curr;
  int index_max;
  IntegerVector Z_max = wrap(Z), Z_last(n),
    phi_max = wrap(phi), phi_last(n),
    K_vec(iter), k_vec(iter);
  NumericVector theta_vec(iter), gamma_vec(iter), alpha_vec(iter),
    lpost_vec(iter), entropy_vec(iter), finite_vec(iter),
    a0_vec(iter), b0_vec(iter);
  running_stat<double> finite_stat, k_stat;
  running_stat_vec<vec> A_stat;
  mat A_mat(n, n, fill::zeros);
  // 04) run
  int g, h, i, j, l, m, o, p, q, s, t;
  // e, f, x, y, z unused
  // A, B, D, J, O, R, S, T, W, X unused
  for (t = 0; t < iter * thin + burn; t++) {
    if (finite) {
      alpha = -gamma;
    }

    // a) update Z
    if (node < 0) {
      stop("gvs: 'node' has to be positive integer.");
    }
    else if (node > n) {
      stop("gvs: 'node' has to be smaller than # of nodes.");
    }
    else if (node == 0) {
      nodes = seqn_rcpp;
    }
    else {
      nodes = Rcpp::RcppArmadillo::sample(seqn_rcpp, node, false, NumericVector::create());
    }
    for (g = 0; g < nodes.size(); g++) {
      q = nodes[g];
      add_group = false;
      p = sigma[q]; // phi[p] == q
      i = Z[p];
      Yrp = Yt.col(p).t();
      Ycp = Y.col(p);
      rp = F.row(p);
      cp = G.col(p);
      rq = U.row(q);
      cq = V.col(q);
      NumericVector w(K);
      IntegerVector seqK = sl1(K);
      P = E + a0;
      P.row(i) -= rp;
      P.col(i) -= cp;
      Q = M + b0;
      Q.row(i) -= rq;
      Q.col(i) -= cq;
      for (j = 0; j < K; j++) {
        w[j] = log(N[j] - alpha) + lincr(P, Q, rp, cp, rq, cq, j);
      }
      w[i] += (log(N[i] - 1.0 - alpha) - log(N[i] - alpha)); // p has to leave i & rejoin
      if (N[i] == 1.0) {
        // Kmp = K - 1
        // symbolically rm C.{col,row}(i) b4 updating
        if (finite) {
          theta = k * (-alpha);
        }
        // adjust weight as if node p joins "new" group i
        w[i] = 0.0
          + (lpg(rp + a0) - dot(rp + a0, log(rq + b0)))
          + (lpg(cp + a0) - dot(cp + a0, log(cq + b0)))
          + log(theta + alpha * (K - 1.0))
          - (lgamma(a0) - a0 * log(b0))
          * (2.0 * K); //// 2022-01-25: works but needs to understand why not 2K-1 nor 2(K-1)
        j = sample_w(seqK, exp(w - max(w)));
        if (j != i) {
          Z[p] = j;
          N[i] -= 1.0;
          N[j] += 1.0;
          E.row(i) -= rp;
          E.col(i) -= cp;
          E.row(j) += rp;
          E.col(j) += cp;
          M.row(i) -= rq;
          M.col(i) -= cq;
          M.row(j) += rq;
          M.col(j) += cq;
          M(i, i) = nC2_weighted(Z, xi, i, dag);
          M(j, j) = nC2_weighted(Z, xi, j, dag);
          F.col(i) -= Ycp;
          F.col(j) += Ycp;
          G.row(i) -= Yrp;
          G.row(j) += Yrp;
          if (dag) {
            if (q != 0) {
              xic = xi_star[q] * xi_star(span(0, q-1));
              U(span(0, q-1), i) -= xic;
              U(span(0, q-1), j) += xic;
            }
            if (q != n-1) {
              xir = xi_star[q] * xi_star(span(q+1, n-1)).t();
              V(i, span(q+1, n-1)) -= xir;
              V(j, span(q+1, n-1)) += xir;
            }
          }
          else {
            Z_star[q] = j; // Z_star needs manual updating for updating xi later
            xic = xi_star[q] * xi_star;
            xic[q] = 0.0;
            U.col(i) -= xic;
            U.col(j) += xic;
            xir = xic.t();
            V.row(i) -= xir;
            V.row(j) += xir;
          }
          // extra lines due to group removal
          Rcout << "Iter " << t + 1 << ", ";
          Rcout << "p = " << p + 1 << ", ";
          Rcout << "- group " << i + 1 << " of " << K << endl;
          rp.shed_col(i);
          cp.shed_row(i);
          rq.shed_col(i);
          cq.shed_row(i);
          N.shed_row(i);
          E.shed_col(i);
          E.shed_row(i);
          M.shed_col(i);
          M.shed_row(i);
          F.shed_col(i);
          G.shed_row(i);
          H.shed_col(i);
          U.shed_col(i);
          I.shed_row(i);
          V.shed_row(i);
          P.shed_col(i);
          P.shed_row(i);
          Q.shed_col(i);
          Q.shed_row(i);
          if (i != K - 1) {
            for (s = i + 1; s < K; s++) {
              u = find(Z == s);
              Z(u).fill(s - 1);
            }
            if (!dag) {
              Z_star = Z(sigma);
            }
          }
          K--;
          lpost_curr += (w[j] - w[i]);
        }
      }
      else {
        if (!finite || (finite && k > K)) {
          // expand weights for potential new group
          if (finite) theta = k * (-alpha);
          ldiff = 0.0 
            + (lpg(rp + a0) - dot(rp + a0, log(rq + b0)))
            + (lpg(cp + a0) - dot(cp + a0, log(cq + b0)))
            + log(theta + alpha * K)
            - (lgamma(a0) - a0 * log(b0))
            * (2.0 * K); //// 2022-01-25: makes sense; c.f. above similar comment
          w.insert(w.end(), ldiff);
          seqK.insert(seqK.end(), K);
        }
        j = sample_w(seqK, exp(w - max(w)));
        if (j != i) { // no update if same group
          if (j == K) { // new group
            add_group = true;
            // extra lines due to new group creation
            Rcout << "Iter " << t + 1 << ", ";
            Rcout << "p = " << p + 1 << ", ";
            Rcout << "+ group " << K + 1 << endl;
            rp.resize(K + 1);
            cp.resize(K + 1);
            rq.resize(K + 1);
            cq.resize(K + 1);
            N.resize(K + 1);
            E.resize(K + 1, K + 1);
            M.resize(K + 1, K + 1);
            F.insert_cols(K, 1);
            G.insert_rows(K, 1);
            H.insert_cols(K, 1);
            U.insert_cols(K, 1);
            I.insert_rows(K, 1);
            V.insert_rows(K, 1);
            P.resize(K + 1, K + 1);
            Q.resize(K + 1, K + 1);
            K++; // correct order - think!
          }
          Z[p] = j;
          N[i] -= 1.0;
          N[j] += 1.0;
          E.row(i) -= rp;
          E.col(i) -= cp;
          E.row(j) += rp;
          E.col(j) += cp;
          M.row(i) -= rq;
          M.col(i) -= cq;
          M.row(j) += rq;
          M.col(j) += cq;
          M(i, i) = nC2_weighted(Z, xi, i, dag);
          M(j, j) = nC2_weighted(Z, xi, j, dag);
          F.col(i) -= Ycp;
          F.col(j) += Ycp;
          G.row(i) -= Yrp;
          G.row(j) += Yrp;
          if (dag) {
            if (q != 0) {
              xic = xi_star[q] * xi_star(span(0, q-1));
              U(span(0, q-1), i) -= xic;
              U(span(0, q-1), j) += xic;
            }
            if (q != n-1) {
              xir = xi_star[q] * xi_star(span(q+1, n-1)).t();
              V(i, span(q+1, n-1)) -= xir;
              V(j, span(q+1, n-1)) += xir;
            }
          }
          else {
            Z_star[q] = j; // Z_star needs manual updating for updating xi later
            xic = xi_star[q] * xi_star;
            xic[q] = 0.0;
            U.col(i) -= xic;
            U.col(j) += xic;
            xir = xic.t();
            V.row(i) -= xir;
            V.row(j) += xir;
          }
          lpost_curr += (w[j] - w[i]);
        }
      }
    }
    H = U.rows(phi);
    I = V.cols(phi);

    /*
    // b) split-merge
    pq = Rcpp::RcppArmadillo::sample(seqn_rcpp, 2, false, NumericVector::create());
    // propose to split into groups i & K
    if (Z[pq[0]] == Z[pq[1]]) { // propose to split
      p = pq[0];
      q = pq[1];
      Ycp = Y.col(p);
      Yrp = Yt.col(p).t();
      i = Z[p];
      node_set = find(Z == i);
      indices = find(node_set != p && node_set != q);
      if (indices.size() > 0 && k > K) {
        // k < K not possible (checks in lpost_par())
        // k == K  =>  k < K + 1  =>  can't split group
        node_set = node_set(indices);
        l = phi[p];
        // set everything equal to before first
        Z0 = Z; // Z_launch
        N0 = N;
        E0 = E;
        M0 = M;
        F0 = F;
        G0 = G;
        U0 = U;
        V0 = V;
        // augment due to p creating new gp
        N0.resize(K + 1);
        E0.resize(K + 1, K + 1);
        M0.resize(K + 1, K + 1);
        F0.insert_cols(K, 1);
        G0.insert_rows(K, 1);
        U0.insert_cols(K, 1);
        V0.insert_rows(K, 1);
        rp = F0.row(p);
        cp = G0.col(p);
        rq = U0.row(l);
        cq = V0.col(l);
        // update due to p moving to new gp
        Z0[p] = K; // 0-indexing
        N0[i] -= 1.0;
        N0[K] += 1.0;
        E0.row(i) -= rp;
        E0.col(i) -= cp;
        E0.row(K) += rp;
        E0.col(K) += cp;
        M0.row(i) -= rq;
        M0.col(i) -= cq;
        M0.row(K) += rq;
        M0.col(K) += cq;
        M0(i, i) = nC2_weighted(Z0, xi, i, dag);
        M0(K, K) = nC2_weighted(Z0, xi, K, dag);
        F0.col(i) -= Ycp;
        F0.col(K) += Ycp;
        G0.row(i) -= Yrp;
        G0.row(K) += Yrp;
        if (l != 0) {
          xic = xi_star[l] * xi_star(span(0, l-1));
          U0(span(0, l-1), i) -= xic;
          U0(span(0, l-1), K) += xic;
        }
        if (l != n-1) {
          xir = xi_star[l] * xi_star(span(l+1, n-1)).t();
          V0(i, span(l+1, n-1)) -= xir;
          V0(K, span(l+1, n-1)) += xir;
        }
        // Gibbs scan to obtain Z1 (Z_split)
        for (h = 0; h < 1 + scan + 1; h++) {
          // 1 + for unbiased scan to obtain Z0 (Z_launch)
          // + 1 for the last scan from Z0 to Z1
          a = 0.0; // log(q(Z1|Z0))
          for (s = 0; s < node_set.size(); s++) {
            o = node_set[s];
            Yco = Y.col(o);
            Yro = Yt.col(o).t();
            g = Z0[o]; // should be all i if h == 0
            l = phi[o];
            rp = F0.row(o);
            cp = G0.col(o);
            rq = U0.row(l);
            cq = V0.col(l);
            if (h == 0) {
              // 1st scan is un{biased/iform/informative}
              lA0 = 0.0;
              lA1 = 0.0;
            }
            else {
              P = E0 + a0;
              P.row(g) -= rp;
              P.col(g) -= cp;
              Q = M0 + b0;
              Q.row(g) -= rq;
              Q.col(g) -= cq;
              lA0 = log(N0[i] - (double) (g == i) - alpha) + lincr(P, Q, rp, cp, rq, cq, i); // group i
              lA1 = log(N0[K] - (double) (g == K) - alpha) + lincr(P, Q, rp, cp, rq, cq, K); // group K(+1)
            }
            r = lr1();
            b = 1.0 / (1.0 + exp(lA1 - lA0));
            if (r < log(b)) {
              // group i
              a += log(b);
              if (g == K) {
                // move from group K(+1) to i
                Z0[o] = i;
                N0[i] += 1.0;
                N0[K] -= 1.0;
                E0.row(i) += rp;
                E0.col(i) += cp;
                E0.row(K) -= rp;
                E0.col(K) -= cp;
                M0.row(i) += rq;
                M0.col(i) += cq;
                M0.row(K) -= rq;
                M0.col(K) -= cq;
                M0(i, i) = nC2_weighted(Z0, xi, i, dag);
                M0(K, K) = nC2_weighted(Z0, xi, K, dag);
                F0.col(i) += Yco;
                F0.col(K) -= Yco;
                G0.row(i) += Yro;
                G0.row(K) -= Yro;
                if (l != 0) {
                  xic = xi_star[l] * xi_star(span(0, l-1));
                  U0(span(0, l-1), i) += xic;
                  U0(span(0, l-1), K) -= xic;
                }
                if (l != n-1) {
                  xir = xi_star[l] * xi_star(span(l+1, n-1)).t();
                  V0(i, span(l+1, n-1)) += xir;
                  V0(K, span(l+1, n-1)) -= xir;
                }
              }
            }
            else {
              // group K(+1)
              a += log(1.0 - b);
              if (g == i) {
                // move from group i to K(+1)
                Z0[o] = K;
                N0[i] -= 1.0;
                N0[K] += 1.0;
                E0.row(i) -= rp;
                E0.col(i) -= cp;
                E0.row(K) += rp;
                E0.col(K) += cp;
                M0.row(i) -= rq;
                M0.col(i) -= cq;
                M0.row(K) += rq;
                M0.col(K) += cq;
                M0(i, i) = nC2_weighted(Z0, xi, i, dag);
                M0(K, K) = nC2_weighted(Z0, xi, K, dag);
                F0.col(i) -= Yco;
                F0.col(K) += Yco;
                G0.row(i) -= Yro;
                G0.row(K) += Yro;
                if (l != 0) {
                  xic = xi_star[l] * xi_star(span(0, l-1));
                  U0(span(0, l-1), i) -= xic;
                  U0(span(0, l-1), K) += xic;
                }
                if (l != n-1) {
                  xir = xi_star[l] * xi_star(span(l+1, n-1)).t();
                  V0(i, span(l+1, n-1)) -= xir;
                  V0(K, span(l+1, n-1)) += xir;
                }
              }
            }
          }
          if (h == scan + 1) {
            Z1 = Z0; // Z_split = Z_launch
            N1 = N0;
            E1 = E0;
            M1 = M0;
            F1 = F0;
            G1 = G0;
            U1 = U0;
            V1 = V0;
          }
        }
        ldiff = lpost0(N1, K + 1, xi, k, gamma, theta, alpha, a0, b0, finite)
          - (lgamma(a0) - a0 * log(b0)) * pow(K + 1.0, 2.0)
          + lpg(E1 + a0) - accu((E1 + a0) % log(M1 + b0))
          - lpost_curr; // computed only once per iteration
        if (lr1() < ldiff - a) {
          Rcout << "Iter " << t + 1 << ", ";
          Rcout << "group " << i + 1 << " -> groups " << i + 1 << " & " << K + 1 << endl;
          Z = Z1;
          N.resize(K + 1);
          E.resize(K + 1, K + 1);
          M.resize(K + 1, K + 1);
          F.insert_cols(K, 1);
          G.insert_rows(K, 1);
          U.insert_cols(K, 1);
          V.insert_rows(K, 1);
          N = N1;
          E = E1;
          M = M1;
          F = F1;
          G = G1;
          U = U1;
          V = V1;
          lpost_curr += ldiff;
          K++;
        }
      }
    }

    // propose to merge groups i & j
    else { // propose to merge groups i & j
      p = pq[0];
      q = pq[1];
      i = Z[p];
      j = Z[q];
      node_set = find(Z == i || Z == j);
      indices = find(node_set != p && node_set != q);
      if (indices.size() > 0) {
        node_set = node_set(indices);
        // set everything equal to before first
        Z0 = Z; // Z_launch
        N0 = N;
        E0 = E;
        M0 = M;
        F0 = F;
        G0 = G;
        U0 = U;
        V0 = V;
        for (h = 0; h < 1 + scan + 1; h++) {
          // 1 + for unbiased scan to obtain Z0 (Z_launch)
          // + 1 for the last scan from Z0 to Z
          a = 0.0; // log(q(Z|Z1)), w/ help of Z0
          for (s = 0; s < node_set.size(); s++) {
            o = node_set[s];
            Yco = Y.col(o);
            Yro = Yt.col(o).t();
            g = Z0[o]; // i or j
            l = phi[o];
            rp = F0.row(o);
            cp = G0.col(o);
            rq = U0.row(l);
            cq = V0.col(l);
            if (h == 0) {
              // 1st scan is un{biased/iform/informative}
              lA0 = 0.0;
              lA1 = 0.0;
            }
            else {
              P = E0 + a0;
              P.row(g) -= rp;
              P.col(g) -= cp;
              Q = M0 + b0;
              Q.row(g) -= rq;
              Q.col(g) -= cq;
              lA0 = log(N0[i] - (double) (g == i) - alpha) + lincr(P, Q, rp, cp, rq, cq, i); // group i
              lA1 = log(N0[j] - (double) (g == j) - alpha) + lincr(P, Q, rp, cp, rq, cq, j); // group j
            }
            r = lr1();
            b = 1.0 / (1.0 + exp(lA1 - lA0));
            if ((h < scan + 1 && r < log(b)) || (h == scan + 1 && Z[o] == i)) {
              // group i
              a += log(b);
              if (g == j) {
                // move from group j to i
                Z0[o] = i;
                N0[i] += 1.0;
                N0[j] -= 1.0;
                E0.row(i) += rp;
                E0.col(i) += cp;
                E0.row(j) -= rp;
                E0.col(j) -= cp;
                M0.row(i) += rq;
                M0.col(i) += cq;
                M0.row(j) -= rq;
                M0.col(j) -= cq;
                M0(i, i) = nC2_weighted(Z0, xi, i, dag);
                M0(j, j) = nC2_weighted(Z0, xi, j, dag);
                F0.col(i) += Yco;
                F0.col(j) -= Yco;
                G0.row(i) += Yro;
                G0.row(j) -= Yro;
                if (l != 0) {
                  xic = xi_star[l] * xi_star(span(0, l-1));
                  U0(span(0, l-1), i) += xic;
                  U0(span(0, l-1), j) -= xic;
                }
                if (l != n-1) {
                  xir = xi_star[l] * xi_star(span(l+1, n-1)).t();
                  V0(i, span(l+1, n-1)) += xir;
                  V0(j, span(l+1, n-1)) -= xir;
                }
              }
            }
            else if ((h < scan + 1 && r >= log(b)) || (h == scan + 1 && Z[o] == j)){
              // group j
              a += log(1.0 - b);
              if (g == i) {
                // move from group i to j
                Z0[o] = j;
                N0[i] -= 1.0;
                N0[j] += 1.0;
                E0.row(i) -= rp;
                E0.col(i) -= cp;
                E0.row(j) += rp;
                E0.col(j) += cp;
                M0.row(i) -= rq;
                M0.col(i) -= cq;
                M0.row(j) += rq;
                M0.col(j) += cq;
                M0(i, i) = nC2_weighted(Z0, xi, i, dag);
                M0(j, j) = nC2_weighted(Z0, xi, j, dag);
                F0.col(i) -= Yco;
                F0.col(j) += Yco;
                G0.row(i) -= Yro;
                G0.row(j) += Yro;
                if (l != 0) {
                  xic = xi_star[l] * xi_star(span(0, l-1));
                  U0(span(0, l-1), i) -= xic;
                  U0(span(0, l-1), j) += xic;
                }
                if (l != n-1) {
                  xir = xi_star[l] * xi_star(span(l+1, n-1)).t();
                  V0(i, span(l+1, n-1)) -= xir;
                  V0(j, span(l+1, n-1)) += xir;
                }
              }
            }
            else {
              Rcout << "Iter " << t + 1 << ", ";
              Rcout << "h = " << h + 1 << ", s = " << s + 1 << endl;
              Rcout << "Attempt to merge but there's escaped case" << endl;
            }
          }
        }
        // set everything equal to before first
        Z1 = Z; // Z_merge
        N1 = N;
        E1 = E;
        M1 = M;
        F1 = F;
        G1 = G;
        U1 = U;
        V1 = V;
        // trim due to all group i nodes moving to j
        node_set = find(Z1 == i); // including p
        for (s = 0; s < node_set.size(); s++) {
          o = node_set[s];
          Yco = Y.col(o);
          Yro = Yt.col(o).t();
          l = phi[o];
          rp = F1.row(o);
          cp = G1.col(o);
          rq = U1.row(l);
          cq = V1.col(l);
          Z1[o] = j;
          N1[i] -= 1.0;
          N1[j] += 1.0;
          E1.row(i) -= rp;
          E1.col(i) -= cp;
          E1.row(j) += rp;
          E1.col(j) += cp;
          M1.row(i) -= rq;
          M1.col(i) -= cq;
          M1.row(j) += rq;
          M1.col(j) += cq;
          M1(i, i) = nC2_weighted(Z1, xi, i, dag);
          M1(j, j) = nC2_weighted(Z1, xi, j, dag);
          F1.col(i) -= Yco;
          F1.col(j) += Yco;
          G1.row(i) -= Yro;
          G1.row(j) += Yro;
          if (l != 0) {
            xic = xi_star[l] * xi_star(span(0, l-1));
            U1(span(0, l-1), i) -= xic;
            U1(span(0, l-1), j) += xic;
          }
          if (l != n-1) {
            xir = xi_star[l] * xi_star(span(l+1, n-1)).t();
            V1(i, span(l+1, n-1)) -= xir;
            V1(j, span(l+1, n-1)) += xir;
          }
        }
        // shed everything after moving
        N1.shed_row(i);
        E1.shed_col(i);
        E1.shed_row(i);
        M1.shed_col(i);
        M1.shed_row(i);
        F1.shed_col(i);
        G1.shed_row(i);
        U1.shed_col(i);
        V1.shed_row(i);
        if (i != K - 1) {
          for (s = i + 1; s < K; s++) {
            u = find(Z1 == s);
            Z1(u).fill(s - 1);
          }
        }
        ldiff = lpost0(N1, K - 1, xi, k, gamma, theta, alpha, a0, b0, finite)
          - (lgamma(a0) - a0 * log(b0)) * pow(K - 1.0, 2.0)
          + lpg(E1 + a0) - accu((E1 + a0) % log(M1 + b0))
          - lpost_curr; // computed only once per iteration
        if (lr1() < ldiff + a) {
          Rcout << "Iter " << t + 1 << ", ";
          Rcout << "groups " << i + 1 << " & " << j + 1 << " -> group " << j + 1 << endl;
          Z = Z1;
          N.shed_row(i);
          E.shed_col(i);
          E.shed_row(i);
          M.shed_col(i);
          M.shed_row(i);
          F.shed_col(i);
          G.shed_row(i);
          U.shed_col(i);
          V.shed_row(i);
          N = N1;
          E = E1;
          M = M1;
          F = F1;
          G = G1;
          U = U1;
          V = V1;
          lpost_curr += ldiff;
          K--;
        }
      }
    }
    H = U.rows(phi);
    I = V.cols(phi);
    */

    // b) update sigma (and/or swap groups?)
    if (dag) {
      for (p = 0; p < n; p++) {
        q = phi[p]; // current position
        m = sample_1(seqL_rcpp); // distance
        g = mod(q + m + n, n); // proposed position
        d = wrap(sigma); // IntegerVector
        d.erase(q);
        d.insert(d.begin() + g, p);
        sigma_prop = as<uvec>(d);
        Y0 = Yt.cols(sigma_prop).t();
        Y0 = Y0.cols(sigma_prop);
        if (Y0.is_trimatu()) {
          // no need to compute if not a topo order
          modulo = false;
          M_prop = M;
          i = Z[p];
          if (q + m < 0 || q + m >= n) {
            // go beyond == append from sigma's other end
            modulo = true;
            m += (q + m < 0) ? n : -n;
          }
          if (m > 0) {
            for (s = q+1; s < q+m+1; s++) {
              o = sigma[s];
              j = Z[o];
              r = xi[p] * xi[o];
              M_prop(i, j) -= r;
              M_prop(j, i) += r;
            }
          }
          else {
            for (s = q-1; s >= q+m; s--) {
              o = sigma[s];
              j = Z[o];
              r = xi[p] * xi[o];
              M_prop(i, j) += r;
              M_prop(j, i) -= r;
            }
          }
          ldiff =
            dot(E.row(i) + a0, (log(M.row(i) + b0) - log(M_prop.row(i) + b0))) +
            dot(E.col(i) + a0, (log(M.col(i) + b0) - log(M_prop.col(i) + b0))) -
            (E(i, i) + a0) * (log(M(i, i) + b0) - log(M_prop(i, i) + b0));
          if (lr1() < ldiff) {
            phi[p] = q + m;
            if (m > 0) {
              for (s = q+1; s < q+m+1; s++) {
                o = sigma[s]; // sigma not updated yet
                j = Z[o];
                r = xi[o] * xi[p];
                H(p, j) -= r;
                H(o, i) += r;
                I(j, p) += r;
                I(i, o) -= r;
                phi[o]--;
              }
            }
            else {
              for (s = q-1; s >= q+m; s--) {
                o = sigma[s]; // sigma not updated yet
                j = Z[o];
                r = xi[o] * xi[p];
                H(p, j) += r;
                H(o, i) -= r;
                I(j, p) -= r;
                I(i, o) += r;
                phi[o]++;
              }
            }
            M = M_prop;
            sigma = sigma_prop;
            lpost_curr += ldiff;
            if (modulo) {
              Rcout << "Iter " << t + 1 << ", ";
              Rcout << "p = " << p + 1 << ", ";
              Rcout << "pos " << q + 1;
              m += (m > 0) ? -n : n; // orig m sampled
              Rcout << ((m > 0) ? " + " : " - ") << abs(m);
              Rcout << " -> " << g + 1 << endl;
            }
          }
        }
      }
      U = H.rows(sigma);
      V = I.cols(sigma);
      xi_star = xi(sigma);
      Z_star = Z(sigma);
    }

    // c) update xi
    for (p = 0; p < n; p++) {
      xip = rnorm(1, xi[p], s_xi[p])[0];
      if (xip > 0.0) {
        xi_prop = xi;
        xi_prop[p] = xip;
        i = Z[p];
        q = phi[p]; // sigma[q] == p
        U0 = U; // prop
        V0 = V; // prop
        if (dag) {
          if (q != 0) {
            U0(span(0, q-1), i) += xi_star(span(0, q-1)) * (xip - xi[p]);
          }
          if (q != n-1) {
            V0(i, span(q+1, n-1)) += xi_star(span(q+1, n-1)).t() * (xip - xi[p]);
          }
        }
        else {
          xic = xi_star * (xip - xi[p]);
          xic[q] = 0.0;
          U0.col(i) += xic;
          xir = xic.t();
          V0.row(i) += xir;
        }
        U0.row(q) *= xip / xi[p];
        V0.col(q) *= xip / xi[p];
        M_prop = M;
        for (j = 0; j < K; j++) {
          u = find(Z_star == j);
          rn = V0.row(i);
          cn = U0.col(i);
          if (j == i && !dag) {
            M_prop(i, i) = nC2_weighted(Z, xi_prop, i, false);
          }
          else {
            M_prop(i, j) = accu(rn(u.t()));
            M_prop(j, i) = accu(cn(u));
          }
        }
        ldiff =
          dot(E.row(i) + a0, (log(M.row(i) + b0) - log(M_prop.row(i) + b0))) +
          dot(E.col(i) + a0, (log(M.col(i) + b0) - log(M_prop.col(i) + b0))) -
          (E(i, i) + a0) * (log(M(i, i) + b0) - log(M_prop(i, i) + b0)) +
          (log(xip) - log(xi[p])) * sY[p] +
          dgamma(tv(xip), axi, 1.0 / bxi, true)[0] -
          dgamma(tv(xi[p]), axi, 1.0 / bxi, true)[0];
        lpost_prop = lpost_curr + ldiff;
        accept_reject = update(xi[p], xip, lpost_curr, lpost_prop, s_xi[p], t, burn);
        if (accept_reject) {
          U = U0;
          V = V0;
          M = M_prop;
          xi_star[q] = xi[p];
        }
      }
    }
    H = U.rows(phi);
    I = V.cols(phi);

    // d) update k / theta
    if (finite) {
      i = rgeom(1, pg)[0] + 1; // 1-based
      i = (lr1() < log(0.5)) ? (k - i) : (k + i);
      ldiff =
        lpost0(N, K, xi, i, gamma, theta, alpha, a0, b0, finite) -
        lpost0(N, K, xi, k, gamma, theta, alpha, a0, b0, finite);
      accept_reject = false;
      if (lr1() < ldiff) {
        accept_reject = true;
        k = i;
        lpost_curr += ldiff;
      }
      k_stat((double) accept_reject);
    }
    else {
      b = theta + rnorm(1, 0.0, s_theta)[0];
      ldiff =
        lpost0(N, K, xi, k, gamma, b, alpha, a0, b0, finite) -
        lpost0(N, K, xi, k, gamma, theta, alpha, a0, b0, finite);
      lpost_prop = lpost_curr + ldiff;
      accept_reject = update(theta, b, lpost_curr, lpost_prop, s_theta, t, burn);
    }
        
    // e) update gamma / alpha
    if (finite) {
      a = gamma * exp(rnorm(1, 0.0, s_gamma)[0]);
      ldiff =
        lpost0(N, K, xi, k, a, theta, alpha, a0, b0, finite) -
        lpost0(N, K, xi, k, gamma, theta, alpha, a0, b0, finite);
      lpost_prop = lpost_curr + ldiff;
      accept_reject = update(gamma, a, lpost_curr, lpost_prop, s_gamma, t, burn, 10.0, log(a) - log(gamma));
      theta = k * gamma;
      alpha = -gamma;
    }
    else { // infinite regime
      a = alpha + rnorm(1, 0.0, s_alpha)[0];
      ldiff =
        lpost0(N, K, xi, k, gamma, theta, a, a0, b0, finite) -
        lpost0(N, K, xi, k, gamma, theta, alpha, a0, b0, finite);
      lpost_prop = lpost_curr + ldiff;
      accept_reject = update(alpha, a, lpost_curr, lpost_prop, s_alpha, t, burn);
    }

    // f) update a0 & b0
    a = a0 + rnorm(1, 0.0, s_a0)[0];
    ldiff =
      lpost(Y, Z, K, sigma, xi, k, gamma, theta, alpha, a, b0, finite) -
      lpost(Y, Z, K, sigma, xi, k, gamma, theta, alpha, a0, b0, finite);
    lpost_prop = lpost_curr + ldiff;
    accept_reject = update(a0, a, lpost_curr, lpost_prop, s_a0, t, burn);
    b = b0 + rnorm(1, 0.0, s_b0)[0];
    ldiff =
      lpost(Y, Z, K, sigma, xi, k, gamma, theta, alpha, a0, b, finite) -
      lpost(Y, Z, K, sigma, xi, k, gamma, theta, alpha, a0, b0, finite);
    lpost_prop = lpost_curr + ldiff;
    accept_reject = update(b0, b, lpost_curr, lpost_prop, s_b0, t, burn);
        
    // g) update regime
    if (p_fin != 0.0 && p_fin != 1.0) {
      if (finite) {
        // sim. theta & alpha from pseudoprior
        a = rgamma(1, a_alpha, 1.0 / b_alpha)[0];
        b = rgamma(1, a_theta, 1.0 / b_theta)[0];
        lpost_fin = lpost0(N, K, xi, k, gamma, theta, alpha, a0, b0, true); // theta & alpha don't matter as finite is true
        lA1 = lpost_fin +
          dgamma(tv(a), a_alpha, 1.0 / b_alpha, true)[0] +
          dgamma(tv(b), a_theta, 1.0 / b_theta, true)[0];
        lpost_inf = lpost0(N, K, xi, k, gamma, b, a, a0, b0, false); // k & gamma don't matter as finite is false
        lA0 = lpost_inf +
          dpois(ti(k), mean_k, true)[0] +
          dgamma(tv(gamma), a_gamma, 1.0 / b_gamma, true)[0];
        if (lr1() > -log(1.0 + exp(lA0 - lA1))) {
          finite = false;
          alpha = a;
          theta = b;
          gamma = -alpha; // just for checking post-run
          lpost_curr += (lpost_inf - lpost_fin);
        }
      }
      else {
        // sim. k & gamma from pseudoprior
        i = rpois(1, mean_k)[0];
        a = rgamma(1, a_gamma, 1.0 / b_gamma)[0];
        lpost_fin = lpost0(N, K, xi, i, a, theta, alpha, a0, b0, true); // theta & alpha don't matter as finite is true
        lA1 = lpost_fin +
          dgamma(tv(alpha), a_alpha, 1.0 / b_alpha, true)[0] +
          dgamma(tv(theta), a_theta, 1.0 / b_theta, true)[0];
        lpost_inf = lpost0(N, K, xi, k, gamma, theta, alpha, a0, b0, false); // k & gamma don't matter as finite is false
        lA0 = lpost_inf +
          dpois(ti(i), mean_k, true)[0] +
          dgamma(tv(a), a_gamma, 1.0 / b_gamma, true)[0];
        if (lr1() < -log(1.0 + exp(lA0 - lA1))) {
          finite = true;
          k = i;
          gamma = a;
          theta = k * gamma; // just for checking post-run
          alpha = -gamma; // just for checking post-run
          lpost_curr += (lpost_fin - lpost_inf);
        }
      }
    }
    finite_stat((double) finite);
    
    // h) print & save
    if ((t + 1) % freq == 0) {
      Rcout << "Iter " << t + 1;
      Rcout << ": Log-posterior = " << lpost_curr << endl;
      Rcout << "K = " << K << endl;
      if (finite) {
        Rcout << "k = " << k << endl;
        Rcout << "k's acceptance rate = " << k_stat.mean() << endl;
        Rcout << "gamma = " << gamma << endl;
        if (t < burn) {
          Rcout << "s_gamma = " << s_gamma << endl;
        }
      }
      else { // infinite regime
        Rcout << "theta = " << theta << endl;
        if (t < burn) {
          Rcout << "s_theta = " << s_theta << endl;
        }
        Rcout << "alpha = " << alpha << endl;
        if (t < burn) {
          Rcout << "s_alpha = " << s_alpha << endl;
        }
      }
      Rcout << "a0 = " << a0 << endl;
      if (t < burn) {
        Rcout << "s_a0 = " << s_a0 << endl;
      }
      Rcout << "b0 = " << b0 << endl;
      if (t < burn) {
        Rcout << "s_b0 = " << s_b0 << endl;
      }
      Rcout << "mean(regime) = " << finite_stat.mean() << endl;
      Rcout << endl;
    }
    if (t >= burn && (t - burn + 1) % thin == 0) {
      s = (t - burn + 1) / thin - 1;
      Z_mat.row(s) = Z.t();
      phi_mat.row(s) = phi.t();
      xi_mat.row(s) = xi.t();
      K_vec[s] = K;
      k_vec[s] = k;
      gamma_vec[s] = gamma;
      theta_vec[s] = theta;
      alpha_vec[s] = alpha;
      a0_vec[s] = a0;
      b0_vec[s] = b0;
      finite_vec[s] = (double) finite;
      entropy_vec[s] = cal_entropy(N, K, n);
      lpost_vec[s] = lpost_curr;
      if (lpost_curr > lpost_max) {
        lpost_max = lpost_curr;
        index_max = s; // 1-indexing later
        Z_max = wrap(Z); // 1-indexing later
        phi_max = wrap(sigma);
      }
      // probs to compare w/ connectivity matrix
      A_mat.zeros();
      for (p = 0; p < n; p++) {
        i = Z[p];
        for (q = 0; q < n; q++) {
          j = Z[q];
          if (q != p) {
            A_mat(p, q) = 1.0 - pow(1.0 - xi[p] * xi[q] / (M(i, j) + b0), E(i, j) - Y(p, q) + a0);
          }
        }
      }
      A_stat(vectorise(A_mat));
    }
  }
  // 05) checks
  Rcout << "Final check: " << endl;
  check_NEM(Y, K, Z, sigma, xi, n, N, E, M, dag);
  lpost_orig = lpost(Y, Z, K, sigma, xi, k, gamma, theta, alpha, a0, b0, finite);
  Rcout << "log-posterior computed from scratch = " << lpost_orig << endl;
  Rcout << "log-posterior updated incrementally = " << lpost_curr << endl;
  Rcout << "difference of the two log-posterior = " << lpost_orig - lpost_curr << endl << endl;
  // 06) return
  // Z, phi, index_max all w/ 1-indexing for R
  scalars["index_max_Z"] = ti(index_max) + 1;
  A_mat = reshape(A_stat.mean(), n, n); // recycle
  DataFrame
    pars =
    DataFrame::create(Named("index") = seq_len(iter),
                      Named("K") = K_vec,
                      Named("k") = k_vec,
                      Named("gamma") = gamma_vec,
                      Named("theta") = theta_vec,
                      Named("alpha") = alpha_vec,
                      Named("a0") = a0_vec,
                      Named("b0") = b0_vec,
                      Named("entropy") = entropy_vec,
                      Named("lpost") = lpost_vec,
                      Named("finite") = finite_vec),
    // all 1-indexing for manipulation in R
    point_est =
    DataFrame::create(Named("Z_last") = Z + 1,
                      Named("phi_last") = phi + 1,
                      Named("xi_last") = xi,
                      Named("Z_max") = Z_max + 1,
                      Named("phi_max") = phi_max + 1,
                      Named("s_xi") = s_xi);
  List output =
    List::create(Named("scalars") = scalars,
                 Named("pars") = pars,
                 Named("point_est") = point_est,
                 Named("Z") = Z_mat + 1,
                 Named("phi") = phi_mat + 1,
                 Named("xi") = xi_mat,
                 Named("A") = A_mat,
                 Named("dag") = dag);
  return output;
}

