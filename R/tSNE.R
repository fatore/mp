#' t-Distributed Stochastic Neighbor Embedding
#'
#' Creates a k-dimensional representation of the data by modeling the
#' probability of picking neighbors using a Gaussian for the high-dimensional
#' data and t-Student for the low-dimensional map and then minimizing the KL
#' divergence between them. This implementation uses the same default parameters
#' as defined on the paper.
#'
#' @param X A data frame, data matrix or dissimilarity matrix.
#' @param Y Initial k-dimensional configuration. If NULL, the method uses a
#' random initial configuration.
#' @param k Target dimension. Avoid anything other than 2 or 3.
#' @param perplexity Perplexity to use (related to the neighborhood size).
#' @param n.iter Number of iterations to perform.
#' @return The kD representation of the data.
#'
#' @references L.J.P. van der Maaten and G.E. Hinton. _Visualizing
#' High-Dimensional Data Using t-SNE._ Journal of Machine Learning Research
#' 9(Nov): 2579-2605, 2008.
#'
#' @examples
#' # Iris example
#' proj = tsne(iris[, 1:4])
#' plot(proj, col = iris$Species)
#'
#' @useDynLib mp
#' @export
tSNE = function(X, Y=NULL, k=2, perplexity=30.0, n.iter=1000) {
  if (!is.matrix(X)) {
    X = as.matrix(X)
  }

  if (is.null(Y)) {
    Y = matrix(runif(nrow(X) * k), ncol = k)
  }

  if (!is.matrix(Y)) {
    Y = as.matrix(Y)
  }

  if (nrow(X) != nrow(Y)) {
    stop("X and Y have different sizes")
  }

  if (ncol(Y) != k) {
    stop("target dimensionality does not match initial map")
  }

  .Call("mp_tSNE", X, Y, perplexity, k, n.iter, PACKAGE="mp")
}
