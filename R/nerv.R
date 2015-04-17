#' Neighbor Retrieval Visualizer (NeRV)
#'
#' Creates a low-dimensional representation of high-dimensional data by
#' optimizing user-defined costs and a tradeoff between mean smoothed precision
#' and mean smoothed recall.
#'
#' @param D A dist object or a symmetric distance matrix.
#' @param lambda The tradeoff between mean smoothed precision and mean smoothed
#'               recall.
#' @return The low-dimensional representation of the data.
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
nerv = function(D, lambda = 0.5) {
  D = as.matrix(D)
  if (!isSymmetric(D)) {
    stop("The distances matrix provided is not symmetric")
  }

  # TODO
}
