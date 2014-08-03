#' Force Scheme
#'
#' Create a 2D representation of the data based on a dissimilarity matrix. A few
#' modifications have been made in relation to the method described in the
#' literature: shuffled indices are used to minimize the order dependency
#' factor, only a fraction of delta is used for better stability and a tolerance
#' factor was introduced as a second stop criterion.
#'
#' @param d A dissimilarity structure such as that returned by dist or a full
#'   symmetric matrix containing the dissimilarities.
#' @param initial A initial 2D configuration (optional). A random configuration
#'   will be created when omitted.
#' @param maxIt The maximum number of iterations that the method will run.
#' @param tol The tolerance for the accumulated error between iterations. Set it
#'   to 0 to guarantee it will run maxIt times.
#' @return The 2D representation of the data.
#'
#' @references Eduardo Tejada, Rosane Minghim, Luis Gustavo Nonato: On improved
#'   projection techniques to support visual exploration of multi-dimensional
#'   data sets. Information Visualization 2(4): 218-231 (2003)
#'
#' @examples
#' result = forceScheme(eurodist)
#' @seealso \code{\link[stats]{dist}} (stats) and \code{\link[proxy]{dist}}
#'   (proxy) for d computation
#'
#' @export
forceScheme = function(d, initial=NULL, maxIt=50, tol=0.1) {
  # define EPSILON (minimum distance value for d2)
  EPSILON = 1E-5

  # define the fraction of delta
  fraction = 8.0

  # convert d to a matrix
  dmat = as.matrix(d)

  # get the number of elements
  n = nrow(dmat)

  # create a random initial configuration if none was given
  if (is.null(initial)) {
    initial = matrix(runif(n*2,0,1), ncol=2)
  }

  p = initial

  # set the previous delta sum as infinity
  prevDeltaSum = 1/0

  # iterate maxIt times
  for (i in 1:maxIt) {
    print(paste("Iteration:", i))

    # create a shuffle array to minimize the order dependency factor
    s_j = sample(n)

    # for each point
    for (j in 1:n) {
      # get the point index
      p1_index = s_j[j]

      # create another shuffle array
      s_k = sample(n)

      # reset delta sum
      deltaSum = 0

      # for each point
      for (k in 1:n) {
        # get the second point index
        p2_index = s_k[k]

        # skip itself
        if (p1_index == p2_index) {
          next;
        }

        # get the distance between the points in R2 (euclidean distance)
        diff = p[p2_index,] - p[p1_index,]
        d2 = sqrt(sum(diff^2))

        # set a minimum distance for d2
        d2 = max(d2, EPSILON)

        # get the distance between the points in RN
        dn = dmat[p1_index, p2_index]

        # calculate delta
        delta = (dn - d2) / fraction

        # move the point_k in the direction of the point_j
        p[p2_index,] = p[p2_index,] + delta * (diff / d2)

        # update delta sum
        deltaSum = deltaSum + abs(delta)
      }
    }

    # check for the tolerance stop criterion
    if (abs(prevDeltaSum - deltaSum) < tol) {
      print(paste("The difference of delta sum is less than tol:", abs(prevDeltaSum - deltaSum)))
      break
    }
    prevDeltaSum = deltaSum
  }

  # return the 2D rojection
  return(p)
}








