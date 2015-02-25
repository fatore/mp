#include <algorithm>
#include <cmath>
#include <RcppArmadillo.h>

static const double ETA = 500;
static const double INITIAL_MOMENTUM   = 0.5;
static const double FINAL_MOMENTUM     = 0.8;
static const double EARLY_EXAGGERATION = 4.;
static const double GAIN_FRACTION      = 0.2;

static const int MOMENTUM_THRESHOLD_ITER     = 20;
static const int EXAGGERATION_THRESHOLD_ITER = 100;
static const int MAX_BINSEARCH_TRIES         = 50;

static void calcP(const arma::mat &X, arma::mat &P, double perplexity, double tol = 1e-5);
static double hBeta(const arma::rowvec &Di, double beta, arma::rowvec &Pi);

double lowerBound(double x) {
    return std::max(x, 1e-12);
}

/*
 * t-SNE C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat tSNE(arma::mat X, arma::mat Y, double perplexity, arma::uword k, arma::uword niter)
{
    arma::uword n = X.n_rows;
    arma::mat P(n, n);
    calcP(X, P, perplexity);
    P = (P + P.t()) / arma::accu(P);
    P *= EARLY_EXAGGERATION;
    P.transform(lowerBound);

    double momentum;

    arma::mat Q(n, n);
    arma::mat dY(n, k),
              gains(n, k),
              iY(n, k);
    for (arma::uword iter = 0; iter < niter; iter++) {
        arma::vec sumY = arma::sum(Y % Y, 1);
        arma::mat num = -2 * (Y * Y.t());
        num.each_col() += sumY;
        arma::inplace_trans(num);
        num.each_col() += sumY;
        num = 1 / (num + 1);
        num.diag() = 0;
        Q = num / arma::accu(num);
        Q.transform(lowerBound); // Q = max(Q, 1e-12);

        for (arma::uword i = 0; i < n; i++) {
            arma::mat tmp = -Y;
            tmp.each_row() += Y.row(i);
            tmp.each_col() %= (P.col(i) - Q.col(i)) % num.col(i);
            dY.row(i) = arma::sum(tmp, 0);
        }

        momentum = (iter < MOMENTUM_THRESHOLD_ITER) ? INITIAL_MOMENTUM : FINAL_MOMENTUM;
        gains = (gains +       GAIN_FRACTION) % ((dY > 0) != (iY > 0))
              + (gains * (1 - GAIN_FRACTION)) % ((dY > 0) == (iY > 0)); // TODO: is this correct?
        iY = momentum * iY - ETA * (gains % dY);
        Y += iY;
        Y.each_row() -= mean(Y, 0);

        if (iter == EXAGGERATION_THRESHOLD_ITER)
            P /= EARLY_EXAGGERATION; // remove early exaggeration
    }

    return Y;
}

static void calcP(const arma::mat &X, arma::mat &P, double perplexity, double tol) {
    arma::colvec sumX = arma::sum(X % X, 1);
    arma::mat D = -2 * (X * X.t());
    D.each_col() += sumX;
    arma::inplace_trans(D);
    D.each_col() += sumX;
    D.diag() = 0;
    double logU = log(perplexity);
    arma::rowvec beta(X.n_rows, arma::fill::ones);

    arma::rowvec Pi(X.n_rows);
    for (arma::uword i = 0; i < X.n_rows; i++) {
        double betaMin = -arma::datum::inf;
        double betaMax =  arma::datum::inf;
        arma::rowvec Di = D.row(i);
        double h = hBeta(Di, beta[i], Pi);

        double hDiff = h - logU;
        for (int tries = 0; fabs(hDiff) > tol && tries < MAX_BINSEARCH_TRIES; tries++) {
            if (hDiff > 0) {
                betaMin = beta[i];
                if (betaMax == arma::datum::inf || betaMax == -arma::datum::inf)
                    beta[i] *= 2;
                else
                    beta[i] = (beta[i] + betaMax) / 2.;
            } else {
                betaMax = beta[i];
                if (betaMin == arma::datum::inf || betaMin == -arma::datum::inf)
                    beta[i] /= 2;
                else
                    beta[i] = (beta[i] + betaMin) / 2.;
            }

            
            h = hBeta(Di, beta[i], Pi);
            hDiff = h - logU;
        }

        P.row(i) = Pi;
    }
}

static inline double hBeta(const arma::rowvec &Di, double beta, arma::rowvec &Pi) {
    Pi = arma::exp(-Di * beta);
    double sumPi = arma::accu(Pi);
    double h = log(sumPi) + beta * arma::accu(Di % Pi) / sumPi;
    Pi /= sumPi;
    return h;
}
