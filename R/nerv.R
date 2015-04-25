#' Neighbor Retrieval Visualizer (NeRV)
#'
#' Creates a k-dimensional representation of data by optimizing a trade-off
#' between mean smoothed precision and mean smoothed recall.
#'
#' @param D A dist object or a symmetric distance matrix.
#' @param lambda The trade-off between mean smoothed precision and mean smoothed
#'               recall.
#' @param Y Initial k-dimensional configuration. Random if omitted.
#' @param k Target dimensionality.
#' @param perplexity A rough upper bound on the neighborhood size.
#' @param max.iter Maximum number of iteration the algorithm performs.
#' @return The k-dimensional representation of the data.
#'
#' @references Venna, J.; Peltonen, J.; Nybo, K.; Aidos, H.; Kaski, S.,
#'   "Information Retrieval Perspective to Nonlinear Dimensionality Reduction for
#'   Data Visualization," Journal of Machine Learning Research, vol.11,
#'   pp.451,490, Feb. 2010.
#'
#' @examples
#'
#' # Iris example
#' emb = nerv(dist(iris[,1:4]))
#' plot(emb, col=iris$Species)
#'
#' @useDynLib mp
#' @export
nerv <- function(D, lambda, Y=NULL, k=2, perplexity=20, max.iter=20) {
  D <- as.matrix(D)
  if (!isSymmetric(D)) {
    stop("The distances matrix provided is not symmetric")
  }

  n <- nrow(D)
  if (is.null(Y)) {
    Y <- matrix(rnorm(n * k), ncol = k)
  }

  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  if (ncol(Y) != k) {
    stop("initial configuration does not match target dimensionality")
  }

  if (nrow(Y) != n) {
    stop("initial configuration does not match dataset size")
  }

  .Call("mp_nerv", D, lambda, Y, k, perplexity, max.iter, Y, PACKAGE="mp")
}
