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

#' Reduce number of stages in MPM simulation output
#'
#' This is designed to make it easy to collapse simulation output from an age-structured
#' model into stages for ease in plotting. It works well with the output provided by pop.projection
#' in the popbio package.
#'
#' @param stage.vectors A matrix, with each row being a class in the MPM and each column being a time
#' step. The stage.vectors element of pop.projection satisfies this.
#' @param stages An optional vector, specifying the stage associated with each class of the MPM. Should
#' be the same length as the number of rows of stage.vectors. If the matrix was produced using a Leslie matrix created using
#'   \code{make_stage4age_matrix}, then this can be left unspecified.
#'
#' @return A matrix similar to the input stage.vectors, but with the rows collapsed into the reduced number of stages.
#' @export
#'
#' @examples
#' # Leslie matrix version of the same data of the loggerhead population
#' A2 <- make_stage4age_matrix(loggerhead[[2]])
#' N0 <- rep(10, nrow(A2))
#' library(popbio)
#' pop.sim <- pop.projection(A2, N0)
#' sv <- collapse_stage.vector(pop.sim$stage.vector)
#' stage.vector.plot(sv)
collapse_stage.vector <- function(stage.vectors, stages = NULL) {
  if (! is.null(stages)){
    stage_names <- stages
  } else {
    stage_names <- unlist(strsplit(rownames(stage.vectors), "[:0-9:]"))
    stage_names <- stage_names[stage_names != ""]
  }
  # Explicitly construct a factor so that the levels aren't sorted alphabetically
  fstage_names <- factor(stage_names, unique(stage_names))
  sv2 <- apply(stage.vectors, 2,
               function(x) as.numeric(by(x, fstage_names, sum, simplify = FALSE)))
  rownames(sv2) <- unique(stage_names)
  sv2
}

dropzeros <- function(A) {
  zerorows <- rowSums(A) == 0
  A[!zerorows, !zerorows]
}