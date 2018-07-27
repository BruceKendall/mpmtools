### eigen-analysis.R
### Asymptotic analysis functions

#' Asymptotic growth rate of an MPM.
#'
#' Calculates lambda_1, the dominant eigenvalue of a matrix, which can be interpreted as
#' the asmyptotic population growth rate of a matrix population model.
#'
#' @param A A square matrix
#'
#' @return The dominant eigenvalue of `A`, which is the asymptotic growth rate of the
#'   deterministic MPM represented by A.
#' @export
#'
#' @examples
#' lambda1(Caswell_Ex_2.1)
lambda1 <- function(A) {
  Re(eigen(A)$values[1])
}