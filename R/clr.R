#' Centred log-ratio (CLR) transformation for compositional data
#'
#' Performs the centred log-ratio (CLR) transform on a numeric matrix or
#' data frame. Optionally handles zeros using \code{zCompositions::cmultRepl}.
#' Rows are treated as observations and columns as parts.
#'
#' For a composition \eqn{X_i = (X_{i1}, \ldots, X_{iD})}, the CLR transform is
#' \deqn{
#'   \mathrm{clr}(X_i)_j = \log(X_{ij}) - \frac{1}{D}\sum_{k=1}^D \log(X_{ik}).
#' }
#' The CLR is scale-invariant, so the rows of \code{X} need not be closed.
#'
#' @param X A numeric matrix or data frame. Rows = observations, columns = parts.
#' @param zero.method Character; how to handle zeros:
#'   \itemize{
#'     \item \code{"none"}: require strictly positive entries; error if any
#'       \eqn{\le 0}.
#'     \item \code{"cmultRepl"}: use \code{zCompositions::cmultRepl} to perform
#'       multiplicative replacement of zeros.
#'   }
#' @param ... Additional arguments passed to \code{zCompositions::cmultRepl}
#'   when \code{zero.method = "cmultRepl"}.
#'
#' @return A numeric matrix of the same dimension as \code{X} containing
#'   CLR-transformed values. Row and column names are preserved if present.
#'
#' @examples
#' ## Example 1: simple positive counts (no zeros)
#' X <- matrix(c(10, 20, 30,
#'               5,  5, 10),
#'             nrow = 2, byrow = TRUE)
#' colnames(X) <- c("A", "B", "C")
#'
#' CLR_X <- clr_transform(X)
#' CLR_X
#'
#' ## Example 2: data with zeros, using zCompositions (if available)
#' X0 <- rbind(c(0, 2, 3),
#'             c(1, 1, 0),
#'             c(1, 1, 5))
#' colnames(X0) <- c("A", "B", "C")
#'
#' CLR_X0 <- clr_transform(X0, zero.method = "cmultRepl")
#' CLR_X0
#'
#' ## Example 3: data frame input also works
#' DF <- data.frame(A = c(4, 1, 2),
#'                  B = c(1, 3, 2),
#'                  C = c(5, 6, 1))
#' CLR_DF <- clr_transform(DF)
#' CLR_DF
#'
#' @export
clr_transform <- function(X,
                          zero.method = c("none", "cmultRepl"),
                          ...) {

  zero.method <- match.arg(zero.method)

  ## --- Input checks ---
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("`X` must be a matrix or data frame.")
  }

  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("`X` must be numeric.")
  }

  if (ncol(X) < 2L) {
    stop("`X` must have at least two columns (parts).")
  }

  if (any(is.na(X))) {
    stop("`X` contains NA values; please handle missing data before CLR.")
  }

  if (any(X < 0)) {
    stop("`X` contains negative values; compositional parts must be >= 0.")
  }

  ## --- Zero handling ---
  if (zero.method == "cmultRepl") {

    if (!requireNamespace("zCompositions", quietly = TRUE)) {
      stop("Package `zCompositions` is required for zero.method = 'cmultRepl'. ",
           "Please install it or use zero.method = 'none'.")
    }

    Xpos <- zCompositions::cmultRepl(X, ...)
    Xpos <- as.matrix(Xpos)

  } else {
    # zero.method == "none"
    if (any(X <= 0)) {
      stop("`X` contains zeros. ",
           "Use zero.method = 'cmultRepl' or replace zeros before CLR.")
    }

    Xpos <- X
  }

  ## --- CLR transform ---
  LogX <- log(Xpos)
  CLR <- LogX - rowMeans(LogX)

  # preserve dimnames explicitly (in case cmultRepl ever messes with them)
  dimnames(CLR) <- dimnames(X)

  return(CLR)
}
