#' Isometric log-ratio (ILR) transformation for compositional data
#'
#' Performs the isometric log-ratio (ILR) transform on a numeric matrix or
#' data frame. Optionally handles zeros using \code{zCompositions::cmultRepl}.
#' Rows are treated as observations and columns as parts.
#'
#' For a composition \eqn{X_i = (X_{i1}, \ldots, X_{iD})}, the ILR transform is
#' defined as an orthonormal rotation of the CLR coordinates:
#' \deqn{
#'   \mathrm{ilr}(X_i) = \mathrm{clr}(X_i)\, \Psi^\top,
#' }
#' where \eqn{\Psi} is a \eqn{(D-1)\times D} contrast matrix with orthonormal
#' rows (i.e., \eqn{\Psi \Psi^\top = I_{D-1}} and each row sums to zero).
#'
#' By default, \code{ilr_transform()} constructs such a \eqn{\Psi} using a scaled
#' Helmert contrast matrix.
#'
#' @param X A numeric matrix or data frame. Rows = observations, columns = parts.
#' @param basis Optional \eqn{(D-1)\times D} numeric matrix giving the ILR
#'   contrast basis \eqn{\Psi}. Rows are ILR coordinates, columns correspond to
#'   parts (columns of \code{X}). If \code{NULL}, a default orthonormal Helmert
#'   basis is constructed.
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
#' @return A numeric matrix with \code{nrow(X)} rows and \code{ncol(X) - 1}
#'   columns containing ILR-transformed values. Row names of \code{X} are
#'   preserved. Column names are taken from \code{rownames(basis)} if present,
#'   otherwise \code{"ilr1"}, \code{"ilr2"}, \dots.
#'
#' @examples
#' ## Example 1: simple positive counts (no zeros)
#' X <- matrix(c(10, 20, 30,
#'               5,  5, 10),
#'             nrow = 2, byrow = TRUE)
#' colnames(X) <- c("A", "B", "C")
#'
#' ILR_X <- ilr_transform(X)
#' ILR_X
#'
#' ## Example 2: data with zeros, using zCompositions (if available)
#' X0 <- rbind(c(0, 2, 3),
#'             c(1, 1, 0),
#'             c(1, 1, 5))
#' colnames(X0) <- c("A", "B", "C")
#'
#' ILR_X0 <- ilr_transform(X0, zero.method = "cmultRepl")
#' ILR_X0
#'
#' ## Example 3: user-supplied ILR basis
#' X3 <- matrix(c(4, 1, 2,
#'                1, 3, 2,
#'                5, 6, 1),
#'              nrow = 3, byrow = TRUE)
#' colnames(X3) <- c("A", "B", "C")
#'
#' ## Simple orthonormal basis for D = 3 (rows sum to zero, rows orthonormal)
#' Psi <- matrix(c( 1, -1,  0,
#'                 1,  1, -2), nrow = 2, byrow = TRUE)
#' Psi <- diag(1 / sqrt(rowSums(Psi^2))) %*% Psi
#' rownames(Psi) <- c("ilr_AB", "ilr_AB_vs_C")
#' colnames(Psi) <- colnames(X3)
#'
#' ILR_X3 <- ilr_transform(X3, basis = Psi)
#' ILR_X3
#'
#' @export
ilr_transform <- function(X,
                          basis = NULL,
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
    stop("`X` contains NA values; please handle missing data before ILR.")
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
           "Use zero.method = 'cmultRepl' or replace zeros before ILR.")
    }

    Xpos <- X
  }

  D <- ncol(Xpos)

  ## --- Build / check ILR basis Ψ ---
  if (is.null(basis)) {
    # Default: scaled Helmert basis (orthonormal rows, row sums = 0)
    # stats::contr.helmert(D) returns a D x (D-1) matrix; transpose to (D-1) x D
    Psi <- t(stats::contr.helmert(D))
    # Normalize rows to unit length
    row_norms <- sqrt(rowSums(Psi^2))
    Psi <- diag(1 / row_norms, nrow = D - 1L) %*% Psi

    # Give some default names
    rownames(Psi) <- paste0("ilr", seq_len(D - 1L))
    colnames(Psi) <- colnames(Xpos)

  } else {
    Psi <- as.matrix(basis)

    if (!is.numeric(Psi)) {
      stop("`basis` must be a numeric matrix.")
    }
    if (!identical(dim(Psi), c(D - 1L, D))) {
      stop("`basis` must have dimension (D-1) x D, where D = ncol(X).")
    }

    # (Optional but cheap) sanity check: rows should sum (approximately) to zero
    row_sums <- rowSums(Psi)
    if (any(abs(row_sums) > 1e-8)) {
      warning("Some rows of `basis` do not sum to zero; ",
              "it may not define a proper ILR contrast.")
    }
  }

  ## --- ILR transform: ilr(X) = clr(X) %*% t(Psi) ---
  LogX <- log(Xpos)
  CLR <- LogX - rowMeans(LogX)

  ILR <- CLR %*% t(Psi)

  # Row names
  rn <- rownames(X)
  if (!is.null(rn) && length(rn) == nrow(ILR)) {
    rownames(ILR) <- rn
  } else {
    rownames(ILR) <- NULL
  }

  # Column names
  bn <- rownames(Psi)
  if (!is.null(bn) && length(bn) == ncol(ILR)) {
    colnames(ILR) <- bn
  } else {
    colnames(ILR) <- paste0("ilr", seq_len(ncol(ILR)))
  }

  return(ILR)
}
