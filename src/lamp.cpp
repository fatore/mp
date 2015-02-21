#include <algorithm>
#include <RcppArmadillo.h>

static const double EPSILON = 1e-3;

/*
 * LAMP C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat lamp(arma::mat X, arma::uvec sampleIndices, arma::mat Ys)
{
    arma::mat Xs = X.rows(sampleIndices);
    arma::uword sampleSize = sampleIndices.n_elem;
    arma::mat projection(X.n_rows, 2);

    for (arma::uword i = 0; i < X.n_rows; i++) {
        arma::rowvec point = X.row(i);

        // calculate alphas
        arma::rowvec alphas(sampleSize);
        for (arma::uword j = 0; j < sampleSize; j++) {
            double dist = arma::accu(arma::square(Xs.row(j) - point));
            alphas[j] = 1. / std::max(dist, EPSILON);
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
        arma::mat M = U.cols(0, 1) * V.t();

        // the projection of point i
        projection.row(i) = (point - Xtil) * M + Ytil;
    }

    return projection;
}
