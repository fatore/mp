#include <RcppArmadillo.h>

/*
 * LSP C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat lsp(const arma::mat & A, const arma::mat & b, const arma::uvec & RsampleIndices, const arma::mat & Ys)
{
    const arma::uvec sampleIndices = RsampleIndices - 1;
    arma::mat projection = arma::solve(A.t() * A, A.t() * b);
    projection.rows(sampleIndices) = Ys;

    return projection;
}
