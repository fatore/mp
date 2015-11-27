#' t-Distributed Stochastic Neighbor Embedding
#'
#' Creates a k-dimensional representation of the data by modeling the
#' probability of picking neighbors using a Gaussian for the high-dimensional
#' data and t-Student for the low-dimensional map and then minimizing the KL
#' divergence between them. This implementation uses the same default parameters
#' as defined by the authors.
#'
#' @param X A data frame, data matrix, dissimilarity (distance) matrix or dist
#'          object.
#' @param Y Initial k-dimensional configuration. If NULL, the method uses a
#'          random initial configuration.
#' @param k Target dimensionality. Avoid anything other than 2 or 3.
#' @param perplexity A rough upper bound on the neighborhood size.
#' @param n.iter Number of iterations to perform.
#' @param eta The "learning rate" for the cost function minimization
#' @param initial.momentum The initial momentum used before changing
#' @param final.momentum The momentum to use on remaining iterations
#' @param early.exaggeration The early exaggeration applied to intial iterations
#' @param gain.fraction
#' @param momentum.threshold.iter Number of iterations before using the final
#'                                momentum
#' @param exaggeration.threshold.iter Number of iterations before using the real
#'                                    probabilities
#' @param max.binsearch.tries Maximum number of tries in binary search for
#'                            parameters to achieve the target perplexity
#'
#' @return The k-dimensional representation of the data.
#'
#' @references L.J.P. van der Maaten and G.E. Hinton. _Visualizing
#' High-Dimensional Data Using t-SNE._ Journal of Machine Learning Research
#' 9(Nov): 2579-2605, 2008.
#'
#' @examples
#' # Iris example
#' emb <- tSNE(iris[, 1:4])
#' plot(emb, col=iris$Species)
#'
#' @useDynLib mp
#' @export
tSNE <- function(X,
                 Y=NULL,
                 k=2,
                 perplexity=30.0,
                 n.iter=1000,
                 eta=500,
                 initial.momentum=0.5,
                 final.momentum=0.8,
                 early.exaggeration=4.0,
                 gain.fraction=0.2,
                 momentum.threshold.iter=20,
                 exaggeration.threshold.iter=100,
                 max.binsearch.tries=50) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  is.dist=F
  tryCatch({
    is.dist <- is.symmetric(X)
  }, error = function(e) {is.dist=F})

  if (is.null(Y)) {
    Y <- matrix(runif(nrow(X) * k), ncol = k)
  }

  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  if (nrow(X) != nrow(Y)) {
    stop("X and Y have different sizes")
  }

  if (ncol(Y) != k) {
    stop("target dimensionality does not match initial map")
  }

  .Call("mp_tSNE",
        X,
        Y,
        perplexity,
        k,
        n.iter,
        is.dist,
        eta,
        initial.momentum,
        final.momentum,
        early.exaggeration,
        gain.fraction,
        momentum.threshold.iter,
        exaggeration.threshold.iter,
        max.binsearch.tries,
        PACKAGE="mp")
}
