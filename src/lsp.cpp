#include <algorithm>
#include <queue>
#include <vector>

#include <RcppArmadillo.h>

static double const EPSILON = 1e-6;

static void dfs(const arma::mat &A, arma::uword i, const std::vector<int> &components);
static void computeNearestNeighbors(const arma::mat &X, arma::mat &A, int k);
static void ensureConnectivity(const arma::mat &X, arma::mat &A);

struct DistComp
{
    bool operator()(const std::pair<arma::uword, double> &x,
                    const std::pair<arma::uword, double> &y) const {
        return x.second < y.second;
    }
};

/*
 * LSP C++ implementation. Refer to the R function for details.
 */
// [[Rcpp::export]]
arma::mat lsp(const arma::mat & X, const arma::uvec & RsampleIndices, const arma::mat & Ys, int k, int q)
{
    const arma::uvec sampleIndices = RsampleIndices - 1;
    arma::uword n = X.n_rows;
    arma::uword nc = sampleIndices.n_elem;
    arma::mat A(n + nc, n, arma::fill::eye);

    computeNearestNeighbors(X, A, k);
    ensureConnectivity(X, A);
    for (arma::uword i = 0; i < nc; i++)
        A(n + i, sampleIndices[i]) = 1;

    arma::mat B(n + nc, q);
    B.submat(0, 0, n - 1, q - 1).zeros();
    B.submat(n, 0, n + nc - 1, q - 1) = Ys;

    arma::mat projection = arma::solve(A.t() * A, A.t() * B);
    return projection;
}

static void computeNearestNeighbors(const arma::mat &X, arma::mat &A, int k)
{
    for (arma::uword i = 0; i < X.n_rows; i++) {
        std::priority_queue<std::pair<arma::uword, double>,
                            std::vector<std::pair<arma::uword, double> >,
                            DistComp> dist;

        const arma::rowvec &ri = X.row(i);
        for (arma::uword j = 0; j < X.n_rows; j++) {
            if (i == j)
                continue;

            double d = arma::norm(ri - X.row(j));
            dist.push(std::make_pair(j, d));
            if (dist.size() > k)
                dist.pop();
        }

        std::vector<std::pair<arma::uword, double> > v;
        for (; !dist.empty(); dist.pop())
            v.push_back(dist.top());

        if (v[0].second < EPSILON)
            A(i, v[0].first) = -1;
        else {
            double sum = 0.;
            for (auto it = v.begin(); it != v.end(); it++) {
                A(i, it->first) = -it->second;
                sum += it->second;
            }

            for (auto it = v.begin(); it != v.end(); it++)
                A(i, it->first) /= sum;
        }
    }
}

static void dfs(const arma::mat &A, arma::uword i, std::vector<int> &components)
{
    for (arma::uword j = 0; j < A.n_cols; j++) {
        if (j == i)
            continue;

        if (components[j] == 0) {
            components[j] = components[i];
            dfs(A, j, components);
        }
    }
}

static void ensureConnectivity(const arma::mat &X, arma::mat &A)
{
    // TODO
    //std::vector<int> components(X.n_rows);

    //int c = 1;
    //for (arma::uword i = 0; i < X.n_rows; i++) {
    //    if (components[i] == 0) {
    //        components[i] = c;
    //        dfs(A, i, components);
    //    }
    //    c++;
    //}
}
