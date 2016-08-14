#include <algorithm>
#include <RcppArmadillo.h>

typedef arma::uvec V;

/*
 * Force Scheme C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat forceScheme(const arma::mat & D,
        arma::mat & Y,
        int max_iter,
        double tol,
        double fraction,
        double EPSILON)
{
    arma::uword n = Y.n_rows;
    V i(n), j(n);
    for (arma::uword k = 0; k < n; k++)
        i[k] = j[k] = k;

    double prev_delta_sum = std::numeric_limits<double>::infinity();
    for (int iter = 0; iter < max_iter; iter++) {
        double delta_sum = 0;

        i = arma::shuffle(i);
        for (V::iterator a = i.begin(); a != i.end(); a++) {
            j = arma::shuffle(j);
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

    return Y;
}
