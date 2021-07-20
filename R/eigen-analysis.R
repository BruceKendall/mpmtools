### eigen-analysis.R
### Asymptotic analysis functions

#' Asymptotic growth rate of an MPM.
#'
#' Calculates lambda_1, the dominant eigenvalue of a matrix, which can be interpreted as
#' the asmyptotic population growth rate of a matrix population model.
#'
#' @param A Either an mpm object or a square matrix representing a matrix
#'   population model
#' @param census A string giving the census to use when extracting the
#'   projection matrix from A (if it's an mpm object). If not specified, the
#'   default census coded into A is used. Value is ignored if A is a matrix.
#'
#' @return The dominant eigenvalue of `A`, which is the asymptotic growth rate of the
#'   deterministic MPM represented by A.
#' @export
#'
#' @examples
#' lambda1(Caswell_Ex_2.1)
lambda1 <- function(A, census = NULL) {
  A <- mpm2A(A, census)
  Re(eigen(A, symmetric = FALSE, only.values = TRUE)$values[1])
}

#' Stable stage distribution of an MPM
#'
#' Calculates w_1, the eigenvector associated with the dominant eigenvector,
#' which can be interpreted as the stable age/stage/size distribution. The
#' vector is scaled to sum to one. Allows stage-for-age Leslie matrices to be
#' reported by stage for ease of visualization.
#'
#' @param A Either an mpm object or a square matrix representing a matrix
#'   population model
#' @param collapse Should classes in the matrix be collapsed into a smaller
#'   number of stage classes?
#' @param stages A vector with length equal to the rank of A giving the stage
#'   name for each matrix class. If A is a Leslie matrix created using
#'   \code{make_stage4age_matrix}, then this can be left unspecified.
#' @param census A string giving the census to use when extracting the
#'   projection matrix from A (if it's an mpm object). If not specified, the
#'   default census coded into A is used. Value is ignored if A is a matrix.
#'
#' @return A vector of the stable stage distribution
#' @export
#'
#' @examples
#' # Reproduce the model in Crowder et al:
#' A1 <- make_stage4age_matrix(loggerhead[[2]], approx_method = "AAS")
#' stable_stage(A1)
#'
#' # Now do the Leslie matrix version of the same data
#' A2 <- make_stage4age_matrix(loggerhead[[2]])
#' stable_stage(A2) # Age structure
#' stable_stage(A2, collapse = TRUE) # Stage structure
stable_stage <- function(A, collapse = FALSE, stages = NULL, census = NULL) {
  A <- mpm2A(A, census)
  w1 <- Re(eigen(A, symmetric = FALSE)$vector)[, 1]
  w1 <- w1 / sum(w1)
  names(w1) <- rownames(A)
  if(collapse) {
    if (! is.null(stages)){
      stage_names <- stages
    } else {
      stage_names <- unlist(strsplit(rownames(A), "[:0-9:]"))
      stage_names <- stage_names[stage_names != ""]
    }
    # Explicitly construct a factor so that the levels aren't sorted alphabetically
    fstage_names <- factor(stage_names, unique(stage_names))

    w1 <- as.numeric(by(w1, fstage_names, sum, simplify = FALSE))
    stage_names <- unique(stage_names)
    names(w1) <- stage_names
  }
  w1
}