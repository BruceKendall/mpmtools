### utils.R
### Utility functions for mpmtools


#' Create Subdiagonal Matrix
#'
#' Create a square matrix with specified subdiagonal elements, or add specified elements
#' to the subdiagonal of an existing matrix.
#'
#' If \code{A} is an integer \eqn{n}, then \code{subdiag(A, sx)} returns an \eqn{n \times
#' n}{n x n} zero matrix with nonzero values on the subdiagonal. If \code{A} is a square
#' matrix, then \code{subdiag(A, sx)} returns the matrix \code{A} with the subdiagonal
#' elements added to it.
#'
#' @param A     Either a positive integer giving the row and column dimensions of the
#'              matrix, or a square matrix
#' @param sx    Values to put on the subdiagonal. Either a scalar (which will be
#'   replicated as needed) or a vector with length one less than `A` or `nrow(A)` (as
#'   appropriate)
#'
#' @return A matrix with specified subdiagonal or with the subdiagonal added to the
#'   input matrix
#' @export
#'
#' @examples
#' # Create zero matrix with non-zero subdiagonal
#' subdiag(5, 1:4)
#'
#' # Add subdiagonal elements to existing matrix
#' A <- matrix(1:16, 4, 4, byrow = TRUE)
#' subdiag(A, 0.5)
#'
#' # Create zero matrix with non-zero subdiagonal that is the size of A
#' subdiag(nrow(A), 0.5)
subdiag <- function(A, sx) {
  if (is.null(dim(A))) {
    A <- matrix(0, A, A)
  }
  stopifnot(is.matrix(A),
            length(sx) == 1 | length(sx) == (nrow(A) - 1),
            nrow(A) == ncol(A))
  n <- nrow(A)
  B <- matrix(0, n, n)
  B[-1, -n] <- diag(sx, n - 1, n - 1)
  return(B + A)
}

dropzeros <- function(A) {
  zerorows <- rowSums(A) == 0
  A[!zerorows, !zerorows]
}