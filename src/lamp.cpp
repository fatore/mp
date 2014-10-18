#include <RcppArmadillo.h>

/*
 * LAMP C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat lamp(arma::mat X, arma::uvec sampleIndices, arma::mat Ys)
{
    Rcpp::Environment svd("package:svd");
    Rcpp::Function propack_svd = svd["propack.svd"];

    arma::mat Xs = X.rows(sampleIndices);

    int sampleSize = sampleIndices.n_elem;
    arma::mat projection(X.n_rows, 2);
    for (int i = 0; i < X.n_rows; i++) {
        arma::rowvec point = X.row(i);

        // calculate alphas
        arma::rowvec alphas(sampleSize);
        for (int j = 0; j < sampleSize; j++)
            alphas[j] = arma::accu(arma::square(Xs.row(j) - point));
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

        // we need only the first 2 vectors in U and V
        Rcpp::List s = Rcpp::as<Rcpp::List>(propack_svd(At * B, Rcpp::Named("neig", 2)));
        arma::mat M = Rcpp::as<arma::mat>(s["u"]) * Rcpp::as<arma::mat>(s["v"]);

        // the projection of point i
        projection.row(i) = (point - Xtil) * M + Ytil;
    }

    return projection;
}
