#' Pekalska's approach to speeding up Sammon's mapping.
#'
#' Creates a k-dimensional representation of the data. As input, a subsample and
#' its k-dimensional mapping are required. The method approximates the subsample
#' mapping to a linear mapping based on the distances matrix of the subsample
#' and then applies the same mapping to all instances.
#'
#' @param D dist object or distances matrix.
#' @param sample.indices The indices of subsamples.
#' @param Ys The subsample mapping (k-dimensional).
#' @return The low-dimensional representation of the data.
#'
#' @references Pekalska, E., de Ridder, D., Duin, R. P., & Kraaijveld, M. A.
#' (1999). A new method of generalizing Sammon mapping with application to
#' algorithm speed-up (pp. 221-228).
#'
#' @useDynLib mp
#' @export
pekalska <- function(D, sample.indices=NULL, Ys=NULL) {
  if (!is.matrix(D)) {
    D <- as.matrix(D)
  }

  n <- nrow(D)

  if (is.null(sample.indices)) {
    sample.indices <- sample(1:n, 3*sqrt(n))
  }

  Ds <- D[sample.indices, sample.indices]

  if (is.null(Ys)) {
    # forceScheme is always 2D
    Ys <- forceScheme(Ds)
  }

  if (!is.matrix(Ys)) {
    Ys <- as.matrix(Ys)
  }

  if (length(sample.indices) != nrow(Ys)) {
    stop("sample.indices and Ys must have the same number of instances")
  }

  Ys <- scale(Ys, center=T, scale=F)
  P <- solve(Ds, Ys)
  Y <- matrix(NA, nrow <- n, ncol <- ncol(Ys))
  Y[sample.indices, ]  <- Ys
  Y[-sample.indices, ] <- D[-sample.indices, sample.indices] %*% P

  Y
}

