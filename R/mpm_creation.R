### mpm_creation.R
### Functions for creating and transforming MPMs

#' Convert prebreeding census MPM to corresponding postbreeding census MPM
#'
#' Prebreeding census MPMs and postbreeding census MPMs are functionally equivalent (they
#' give the same results when analyzed), as the only difference is when the population is
#' counted. It is far easier to construct a prebreeding census MPM without inadvertently
#' introducing errors. However, such a model lacks newborns in the population vector;
#' users may want to include the newborns for ease of comparing the model projections with
#' actual censuses. This function automates the conversion of a prebreeding into a
#' postbreeding census MPM.
#'
#' S0 must be provided. In addition, provide either Amat alone (in which case the top row
#' of Amat will be assumed to comprise fertility coefficients), Amat and Fmat, or Fmat and
#' Umat.
#'
#' @section To do:
#' Account for multiple non-zero rows of Fmat
#'
#' Allow one newborn type to grow to multiple classes at age one
#'
#' Add argument testing
#'
#' Allow Amat to be an mpm object
#'
#' @param S0 Newborn survival. If Fmat is not supplied, or has only one non-zero row, then
#'   S0 should be a scalar. If Fmat has more than one non-zero row, then S0 may be a
#'   scalar (in which case all newborn classes have the same survival) or a vector with
#'   length equalling the number of non-zero rows in Fmat.
#' @param Amat A square matrix representing a prebreeding census MPM
#' @param Fmat A matrix the same size of Amat, containing positive elements for the
#'   fertility coefficients, and zeros elsewhere
#' @param Umat A matrix the same size of Amat, containing positive elements for the
#'   survival/transition coefficients, and zeros elsewhere
#'
#' @return A matrix containing the postbreeding census MPM that corresponds to Amat (or to
#'   Fmat+Umat). Its rank will be one larger than that of Amat
#' @export
#'
#' @examples
pre_to_post <- function(S0, Amat = NULL, Fmat = NULL, Umat = NULL) {
  if (is.null(Fmat)) {
    Fmat <- 0 * Amat
    Fmat[1, ] <- Amat[1, ]
    warning("Fmat not supplied; assuming top row of Amat are fecundity coefficients")
  }

  if (is.null(Umat)) {
    Umat <- Amat - Fmat
  }

  rank_A <- nrow(Amat)
  class1 <- (1:rank_A)[rowSums(Fmat) > 0] # The "yearling" class
  bx <- Fmat[class1, ] / S0

  # Create the ouput U matrix. Currently assumes only one newborn type
  U_post <- Umat
  U_post <- cbind(numeric(rank_A), U_post)
  U_post <- rbind(numeric(rank_A + 1), U_post)
  U_post[class1 + 1, 1] <- S0

  # Create the output F matrix
  ### THIS IS NOT YET CORRECT!!! ###
  F_post <- 0 * U_post
  F_post[1, 1:rank_A] <- bx * rowSums(Umat)

  return(F_post + U_post)
}