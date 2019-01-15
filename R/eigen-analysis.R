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
  Re(eigen(A, symmetric = FALSE, only.values = TRUE)$values[1])
}

stable_stage <- function(A, collapse = FALSE, stages = NULL) {
  w1 <- Re(eigen(A, symmetric = FALSE)$vector)[, 1]
  if(collapse) {
    if (! is.null(stages)){
      stage_names <- stages
    } else {
      stage_names <- unlist(strsplit(rownames(A), "[:0-9:]"))
      stage_names <- stage_names[stage_names != ""]
      stage_names <- unique(stage_names)
    }
    # Aggregate somhow. Maybe using something like by(w1, stage_names, sum) [with the second version of stage_names]
  }
}