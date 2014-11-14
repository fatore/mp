#include <algorithm>
#include <vector>
#include <RcppArmadillo.h>

/*
 * Arrange the N elements of ARRAY in random order.
 * Only effective if N is much smaller than RAND_MAX;
 * if this may not be the case, use a better random
 * number generator.
 */
static void shuffle(std::vector<size_t> &array, size_t n)
{
    if (n > 1) {
        for (size_t i = 0; i < n - 1; i++) {
            size_t j = i + unif_rand() / (1 / (n - i) + 1);
            std::swap(array[i], array[j]);
        }
    }
}

/*
 * Force Scheme C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat forceScheme(arma::mat p,
        arma::mat dmat,
        int max_it,
        double tol,
        double EPSILON,
        double fraction)
{
    // get R random number generator
    GetRNGstate();

    arma::rowvec diff(2);
    size_t size = (size_t) p.n_rows;

    std::vector<size_t> si(size), sj(size);
    for (size_t i = 0; i < size; i++)
        si[i] = sj[i] = i;

    double prev_delta_sum = 1. / 0.;
    for (int iter = 0; iter < max_it; iter++) {
        double delta_sum = 0;
        shuffle(si, size);
        for (size_t i = 0; i < size; i++) {
            size_t p1_index = si[i];
            shuffle(sj, size);
            for (size_t j = 0; j < size; j++) {
                size_t p2_index = sj[j];
                if (p1_index == p2_index)
                    continue;

                diff = p.row(p2_index) - p.row(p1_index);
                // TODO: Verify alternatives to avoid over/underflow
                double d2 = sqrt(arma::accu(diff % diff));
                if (d2 < EPSILON)
                    d2 = EPSILON;

                double delta = (dmat(p1_index, p2_index) - d2) / fraction;
                delta_sum += fabs(delta);
                p.row(p2_index) += delta * (diff / d2);
            }
        }

        if (fabs(prev_delta_sum - delta_sum) < tol)
            break;
        prev_delta_sum = delta_sum;
    }

    // free R random number generator
    PutRNGstate();
    return p;
}
