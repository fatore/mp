#include <algorithm>
#include <RcppArmadillo.h>

static const double EPSILON = 1e-7;

/*
 * PLMP C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat plmp(const arma::mat & X, const arma::uvec & RsampleIndices, const arma::mat & Ys)
{
    // R indices start from 1; we start from 0
    const arma::uvec sampleIndices = RsampleIndices - 1;
    arma::uword sampleSize = sampleIndices.n_elem;
    const arma::mat &Xs_ = X.rows(sampleIndices);

    // Used to avoid singular linear systems
    arma::mat Xs = Xs_.cols(arma::any(Xs_ > EPSILON));
    Xs.diag().ones();

    arma::mat P = arma::solve(Xs.t() * Xs, Xs.t() * Ys);
    arma::mat projection = X * P;
    projection.rows(sampleIndices) = Ys;

    return projection;
}
