#include <algorithm>
#include <cmath>
#include <RcppArmadillo.h>

static const double MIN_GAIN = 1e-2;
static const double EPSILON  = 1e-12;

static void calcP(const arma::mat &X, arma::mat &P, double perplexity, int MAX_BINSEARCH_TRIES, double tol = 1e-5);
static double hBeta(const arma::rowvec &Di, double beta, arma::rowvec &Pi);

class MaxTransform
{
private:
    double m_max;

public:
    MaxTransform(double max): m_max(max) {}

    double operator()(double x) {
        return std::max(x, m_max);
    }

    void setMax(double max) { m_max = max; }
};

/*
 * t-SNE C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat tSNE(const arma::mat & X,
               const arma::mat & initialY,
               double perplexity,
               arma::uword k,
               arma::uword niter,
               bool isDist,
               double ETA,
               double INITIAL_MOMENTUM,
               double FINAL_MOMENTUM,
               double EARLY_EXAGGERATION,
               double GAIN_FRACTION,
               int MOMENTUM_THRESHOLD_ITER,
               int EXAGGERATION_THRESHOLD_ITER,
               int MAX_BINSEARCH_TRIES)
{
    double momentum;
    arma::uword n = X.n_rows;
    arma::mat Q(n, n);
    arma::mat dY(n, k),
              gains(n, k, arma::fill::ones),
              iY(n, k, arma::fill::zeros);
    MaxTransform maxTransform(0);
    arma::mat Y = initialY;

    arma::mat P(n, n, arma::fill::zeros);
    if (!isDist) {
        arma::mat D = -2 * (X * X.t());
        arma::colvec sumX = arma::sum(X % X, 1);
        D.each_col() += sumX;
        arma::inplace_trans(D);
        D.each_col() += sumX;
        D.diag() *= 0;

        calcP(D, P, perplexity, MAX_BINSEARCH_TRIES);
    } else
        calcP(X, P, perplexity, MAX_BINSEARCH_TRIES);
    P = (P + P.t());
    P /= arma::accu(P);
    P *= EARLY_EXAGGERATION;
    maxTransform.setMax(EPSILON);
    P.transform(maxTransform); // P = max(P, EPSILON)

    for (arma::uword iter = 0; iter < niter; iter++) {
        arma::vec sumY = arma::sum(Y % Y, 1);
        arma::mat num = -2. * (Y * Y.t());
        num.each_col() += sumY;
        arma::inplace_trans(num);
        num.each_col() += sumY;
        num = 1. / (1. + num);
        num.diag() *= 0;
        Q = num / arma::accu(num);
        maxTransform.setMax(EPSILON);
        Q.transform(maxTransform); // Q = max(Q, EPSILON);

        for (arma::uword i = 0; i < n; i++) {
            arma::mat tmp = -Y;
            tmp.each_row() += Y.row(i);
            tmp.each_col() %= (P.col(i) - Q.col(i)) % num.col(i);
            dY.row(i) = arma::sum(tmp, 0);
        }

        momentum = (iter < MOMENTUM_THRESHOLD_ITER) ? INITIAL_MOMENTUM : FINAL_MOMENTUM;
        gains = (gains +       GAIN_FRACTION) % ((dY > 0) != (iY > 0))
              + (gains * (1 - GAIN_FRACTION)) % ((dY > 0) == (iY > 0));
        maxTransform.setMax(MIN_GAIN);
        gains.transform(maxTransform); // gains = max(gains, MIN_GAIN)
        iY = momentum * iY - ETA * (gains % dY);
        Y += iY;
        Y.each_row() -= mean(Y, 0);

        if (iter == EXAGGERATION_THRESHOLD_ITER)
            P /= EARLY_EXAGGERATION; // remove early exaggeration
    }

    return Y;
}

static void calcP(const arma::mat &D, arma::mat &P, double perplexity, int MAX_BINSEARCH_TRIES, double tol) {
    double logU = log(perplexity);
    arma::rowvec beta(D.n_rows, arma::fill::ones);

    arma::rowvec Pi(D.n_rows);
    for (arma::uword i = 0; i < D.n_rows; i++) {
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
