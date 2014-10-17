#' Force Scheme Projection
#'
#' Creates a 2D representation of the data based on a dissimilarity matrix. A few
#' modifications have been made in relation to the method described in the
#' literature: shuffled indices are used to minimize the order dependency
#' factor, only a fraction of delta is used for better stability and a tolerance
#' factor was introduced as a second stop criterion.
#'
#' @param d A dissimilarity structure such as that returned by dist or a full
#'   symmetric matrix containing the dissimilarities.
#' @param initial A initial 2D configuration (optional). A random configuration
#'   will be created when omitted.
#' @param max.it The maximum number of iterations that the method will run.
#' @param tol The tolerance for the accumulated error between iterations. Set it
#'   to 0 to guarantee it will run max.it times.
#' @param verbose A flag that indicates if the progress should be printed.
#' @return The 2D representation of the data.
#'
#' @references Eduardo Tejada, Rosane Minghim, Luis Gustavo Nonato: On improved
#'   projection techniques to support visual exploration of multi-dimensional
#'   data sets. Information Visualization 2(4): 218-231 (2003)
#'
#' @examples
#' # Eurodist example
#' proj = forceScheme(eurodist)
#' plot(proj, type = "n", xlab ="", ylab ="", asp=1, axes=FALSE, main="")
#' text(proj, labels(eurodist), cex = 0.6)
#'
#' # Iris example
#' proj = forceScheme(dist(iris[,1:4]))
#' plot(proj, col=iris$Species)
#' @seealso \code{\link[stats]{dist}} (stats) and \code{\link[proxy]{dist}}
#'   (proxy) for d computation
#'
#' @useDynLib mp
#' @export
forceScheme = function(d, initial=NULL, max.it=50, tol=0.1, verbose=F) {
  EPSILON = 1E-5 # minimum distance between points
  fraction = 8.0 # fraction of delta

  # convert d to a matrix
  dmat = as.matrix(d)

  # get the number of elements
  n = nrow(dmat)

  # create a random initial configuration if none was given
  if (is.null(initial)) {
    initial = matrix(runif(n*2,0,1), ncol=2)
  }

  # check initial configuration
  if (nrow(initial) != n) {
    stop("Provided initial configuration does not match the input matrix length");
  }

  p = initial

  # switch core implementation
  pC = .C("force_scheme",
          p = as.numeric(p),
          dmat = as.numeric(dmat),
          n = as.integer(n),
          max_it = as.integer(max.it),
          tol = as.double(tol),
          EPSILON = as.double(EPSILON),
          fraction = as.double(fraction))$p

  p = cbind(pC[1:n], pC[(n+1):(2*n)])

  # return the 2D rojection
  return(p)
}








