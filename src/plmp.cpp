#include <algorithm>
#include <RcppArmadillo.h>

/*
 * PLMP C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat plmp(const arma::mat & X, const arma::uvec & RsampleIndices, const arma::mat & Ys)
{
    // R indices start from 1; we start from 0
    const arma::uvec sampleIndices = RsampleIndices - 1;
    const arma::mat &Xs = X.rows(sampleIndices);
    arma::uword sampleSize = sampleIndices.n_elem;
    arma::mat projection(X.n_rows, Ys.n_cols);

    arma::mat P = arma::solve(Xs.t() * Xs, Xs.t() * Ys);
    projection = X * P;
    projection.rows(sampleIndices) = Ys;

    return projection;
}
