#' Part-Linear Multidimensional Projection
#'
#' Creates a k-dimensional representation of the data. As input, a subsample and
#' its k-dimensional mapping (control points) are required. The method
#' approximates the subsample mapping to a linear mapping and then applies the
#' same mapping to all instances.
#'
#' @param X A dataframe or matrix representing the data.
#' @param sample.indices The indices of subsamples used as control points.
#' @param Ys The control points.
#' @return The low-dimensional representation of the data.
#'
#' @references Paulovich, F.V.; Silva, C.T.; Nonato, L.G., "Two-Phase Mapping
#'   for Projecting Massive Data Sets," Visualization and Computer Graphics,
#'   IEEE Transactions on , vol.16, no.6, pp.1281,1290, Nov.-Dec. 2010.
#'
#' @examples
#'
#' # Iris example
#' emb <- plmp(iris[,1:4])
#' plot(emb, col=iris$Species)
#'
#' @useDynLib mp
#' @export
plmp <- function(X, sample.indices=NULL, Ys=NULL) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  n <- nrow(X)

  if (is.null(sample.indices)) {
    sample.indices <- sample(1:n, 3*sqrt(n))
  }

  if (is.null(Ys)) {
    Ys <- forceScheme(dist(X[sample.indices, ]))
  }

  if (!is.matrix(Ys)) {
    Ys <- as.matrix(Ys)
  }

  if (length(sample.indices) != nrow(Ys)) {
    stop("sample.indices and Ys must have the same number of instances")
  }

  Ys <- scale(Ys, center=T, scale=F)
  .Call("mp_plmp", X, sample.indices, Ys, package="mp")
}

