#include <algorithm>
#include <vector>
#include <RcppArmadillo.h>

typedef std::vector<size_t> V;

/*
 * Fisher-Yates (aka Knuth) shuffle.
 */
static void shuffle(V &array, size_t n)
{
    if (n < 2)
        return;

    for (size_t i = n - 1; i > 0; i--) {
        size_t j = (size_t) (unif_rand() * (i + 1));
        std::swap(array[i], array[j]);
    }
}

/*
 * Force Scheme C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat forceScheme(arma::mat D,
        arma::mat Y,
        int max_iter,
        double tol,
        double fraction,
        double EPSILON)
{
    // get R random number generator
    GetRNGstate();

    size_t n = (size_t) Y.n_rows;
    V i(n), j(n);
    for (size_t k = 0; k < n; k++)
        i[k] = j[k] = k;

    double prev_delta_sum = 1. / 0.;
    for (int iter = 0; iter < max_iter; iter++) {
        double delta_sum = 0;

        shuffle(i, n);
        for (V::iterator a = i.begin(); a != i.end(); a++) {
            shuffle(j, n);
            for (V::iterator b = j.begin(); b != j.end(); b++) {
                if (*a == *b)
                    continue;

                arma::rowvec direction(Y.row(*b) - Y.row(*a));
                double d2 = std::max(arma::norm(direction, 2), EPSILON);
                double delta = (D(*a, *b) - d2) / fraction;
                delta_sum += fabs(delta);
                Y.row(*b) += delta * (direction / d2);
            }
        }

        if (fabs(prev_delta_sum - delta_sum) < tol)
            break;
        prev_delta_sum = delta_sum;
    }

    // free R random number generator
    PutRNGstate();
    return Y;
}
