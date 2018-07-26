### mpm_creation.R
### Functions for creating and transforming MPMs

#' Create a Leslie Matrix Model
#'
#' Create a Leslie matrix from a schedule of age-specific survival and birth rates.
#' The age classes must be evenly spaced, and survival is from one age class to the next.
#' The final age class may be terminal (no further survival) or can be constructed with
#' a self-loop for indefinite age classes having the same demography as the final one.
#'
#' If the survival of the final age class is non-zero, then a self-loop will be added that
#' creates a final indefinite "stage". Note that if the final survival is zero and a
#' postbreeding census model created, then the final column will be all zero. That's ok.
#'
#' @param x Either (1) a data frame containing columns named "x" (age, starting with age
#'   zero), "sx" (containing one-timestep survival rates for each age) and "mx"
#'   (containing birth rates associated with each age); or (2) an integer giving the
#'   maximum age (in timesteps) for the model; or (3) a vector of ages.
#' @param sx A vector of age-specific survival rates. Not needed if x is a data frame
#' @param mx A vector of age-specific birth rates Not needed if x is a data frame
#' @param model The type of model. Currently supported options are "pre" (prebreeding
#'   census model, the default) and "post" (postbreeding census model)
#'
#' @return A matrix containing the Leslie matrix model
#' @export
#'
#' @examples
#' # Create a demography schedule, with juvenile and senescent age classes
#' demog_sched <- data.frame(x = 0:7,
#'                           sx = c(0.05, 0.2, 0.35, 0.8, 0.9, 0.9, 0.75, 0.4),
#'                           mx = c(0, 0, 0, 0.5, 1, 3, 3, 1.5))
#'
#' # Make a Leslie matrix using the data frame and the default model structures
#' make_Leslie_matrix(demog_sched)
#'
#' # Supply the vectors directly, and get a postbreeding census model
#' with(demog_sched, make_Leslie_matrix(0:7, sx, mx, model = "post"))
#'
#' # Supply x as an integer, and get a terminal model ("true" Leslie matrix)
#' with(demog_sched, make_Leslie_matrix(8, c(sx, 0), c(mx, 1)))
make_Leslie_matrix <- function(x, sx = NULL, mx = NULL, model = c("pre", "post")) {
  if(is.data.frame(x)) {
    stopifnot(c("x", "sx", "mx") %in% names(x))
    sx <- x$sx
    mx <- x$mx
    x <- x$x
  }
  if (length(x) == 1) x <- 0:x

  # Make prebreeding census model
  n <- length(x) - 1
  A <- subdiag(n, sx[2:n])
  A[n, n] <- sx[n + 1]
  A[1, ] <- sx[1] * mx[-1]

  if (model[1] == "post") A <- suppressWarnings(pre_to_post(sx[1], A))

  return(A)
}

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
#' Currently, the only matrix structure supported is one in which only the first row of
#' the matrix contains fertility coefficients. This is particularly suitable for
#' non-spatial, one-sex or asexual models of animals (or plants without a seedbank).
#'
#' @section To do:
#' \itemize{
#'   \item Account for multiple non-zero rows of Fmat
#'   \item Allow one newborn type to grow to multiple classes at age one
#'   \item Allow multiple newborn types
#'   \item Add argument testing
#'   \item Allow Amat to be an mpm object
#' }
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

  #-------------------------------------------------------
  # Case 1 (the easy one): only top row of Fmat is nonzero

  if (length(class1) == 1 & class1[1] == 1) {
    # Create the ouput U matrix.
    U_post <- Umat
    U_post <- cbind(0, U_post)
    U_post <- rbind(0, U_post)
    U_post[class1 + 1, 1] <- S0

    # Create the output F matrix
    ### THIS IS NOT YET CORRECT!!! ###
    F_post <- 0 * U_post
    # There is probably a cleverer way to do this...
    F_post[1, 1] <- S0 * bx[1]
    for (i in 2:(rank_A + 1)) {
      F_post[1, i] <- sum(bx * Umat[, i - 1])
    }
  } else {
    stop("The structure of the supplied Fmat is not yet supported")
  }
  return(F_post + U_post)
}