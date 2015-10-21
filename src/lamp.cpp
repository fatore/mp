#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppArmadillo.h>


static const double EPSILON = 1e-6;

/*
 * LAMP C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
arma::mat lamp(const arma::mat & X, const arma::uvec & RsampleIndices, const arma::mat & Ys, double cp)
{
    // R indices start from 1; we start from 0
    const arma::uvec sampleIndices = RsampleIndices - 1;
    const arma::mat &Xs = X.rows(sampleIndices);
    arma::uword sampleSize = sampleIndices.n_elem;
    arma::mat projection(X.n_rows, Ys.n_cols);

    #pragma omp parallel for shared(X, Xs, sampleIndices, Ys, projection)
    for (arma::uword i = 0; i < X.n_rows; i++) {
        const arma::rowvec &point = X.row(i);

        // calculate alphas
        bool skip = false;
        arma::rowvec alphas(sampleSize);
        for (arma::uword j = 0; j < sampleSize; j++) {
            double dist = arma::accu(arma::square(Xs.row(j) - point));
            if (dist < EPSILON) {
                // point is too close to sample point; position them equally
                projection.row(i) = Ys.row(j);
                skip = true;
                break;
            }

            alphas[j] = 1.0 / dist;
        }

        if (skip)
            continue;

        arma::uword c = (arma::uword) (sampleSize * cp);
        if (c < sampleSize) {
            arma::uvec idx = arma::sort_index(alphas, "descend");
            for (arma::uword j = c; j < sampleSize; j++)
                alphas[idx[j]] = 0;
        }

        double alphas_sum = arma::accu(alphas);
        arma::rowvec alphas_sqrt = arma::sqrt(alphas);

        // calculate \tilde{X} and \tilde{Y}
        arma::rowvec Xtil = arma::sum(alphas * Xs, 0) / alphas_sum;
        arma::rowvec Ytil = arma::sum(alphas * Ys, 0) / alphas_sum;

        // calculate \hat{X} and \hat{Y}
        arma::mat Xhat = Xs;
        Xhat.each_row() -= Xtil;
        arma::mat Yhat = Ys;
        Yhat.each_row() -= Ytil;

        // calculate A and B
        arma::mat At = Xhat.t();
        At.each_row() %= alphas_sqrt;
        arma::mat B = Yhat;
        B.each_col() %= alphas_sqrt.t();

        arma::mat U, V;
        arma::vec s;
        arma::svd(U, s, V, At * B);
        arma::mat M = U.head_cols(Ys.n_cols) * V.t();

        projection.row(i) = (point - Xtil) * M + Ytil;
    }

    return projection;
}
