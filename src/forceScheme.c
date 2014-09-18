#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/*
 * Arrange the N elements of ARRAY in random order.
 * Only effective if N is much smaller than RAND_MAX;
 * if this may not be the case, use a better random
 * number generator.
 */
static void shuffle(size_t *array, size_t n)
{
    if (n > 1) {
        size_t i;
        for (i = 0; i < n - 1; i++) {
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}

void force_scheme(double p[],
        const double dmat[],
        const int *n,
        const int *max_it,
        const double *tol,
        const double *EPSILON,
        const double *fraction)
{
    double prev_delta_sum = 1. / 0., delta_sum, d2, dn, delta, diff_x, diff_y;
    size_t i, j, k, p1_index, p2_index;
    size_t size = (size_t) *n;

    size_t s_j[size], s_k[size];
    for (i = 0; i < size; i++)
        s_j[i] = s_k[i] = i;

    int asd = 0;

    for (i = 0; i < *max_it; i++) {
        shuffle(s_j, size);
        for (j = 0; j < size; j++) {
            p1_index = s_j[j];

            delta_sum = 0;
            shuffle(s_k, size);
            for (k = 0; k < size; k++) {
                p2_index = s_k[k];

                if (p1_index == p2_index)
                    continue;

                diff_x = p[p2_index] - p[p1_index];
                diff_y = p[p2_index + size] - p[p1_index + size];

                // TODO: Verify alternatives to avoid over/underflow
                d2 = sqrt(diff_x * diff_x + diff_y * diff_y);

                if (d2 < *EPSILON)
                    d2 = *EPSILON;

                dn = dmat[p1_index * size + p2_index];
                delta = (dn - d2) / *fraction;
                delta_sum += fabs(delta);
                p[p2_index]        += delta * (diff_x / d2);
                p[p2_index + size] += delta * (diff_y / d2);
            }

        }

        if (fabs(prev_delta_sum - delta_sum) < *tol)
            break;
        prev_delta_sum = delta_sum;
    }
}
