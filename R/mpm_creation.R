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

#' Create a Lefkovitch stage-structured matrix
#'
#' In a Lefkovitch model, individuals may either stay in the current stage or transition
#' to the next stage, with the rate of moving from stage i to stage i+1, conditional on
#' survival, being given as g_i (commonly called "growth" but really representing
#' "maturation"). This function creates a Lefkovitch matrix using user-supplied
#' information on maturation rates as well as stage-specific survival and maternity
#' rates.
#'
#' @param stage_table Either a vector of stage names or a data frame with columns
#' "stage_name", "survival", "maternity", and "maturation". In the latter case the next
#' three arguments do not need to be provided
#' @param survival A vector of stage-specific survival, on a per-timestep basis
#' @param maternity A vector of stage-specific maternities (number of offspring produced
#' by an individual in a given stage)
#' @param maturation A vector of the maturation rate, g_i, representing the fraction of
#' individuals that (conditional on survival) should move on to the next stage. If the
#' stage only lasts one timestep (all surviving individuals mature) then use a value
#' of 1.
#' @param model Whether the matrix should represent a prebreeding ("pre") or postbreeding
#' ("post") census model. Defaults to "post"
#'
#' @return A matrix representing the Lefkovitch MPM.
#' @export
#'
#' @examples
make_Lefkovitch_matrix <- function(stage_table, survival = stage_table$survival,
                                   maternity = stage_table$maternity,
                                   maturation = stage_table$maturation,
                                   model = c("post", "pre")) {

  # Argument unpacking
  if(is.data.frame(stage_table)) {
    stage_name <- stage_table$stage
  } else {
    stage_name <- stage_table
  }
  model <- match.arg(model)

  # Useful indicies
  n_stage <- length(stage_name)
  stage_index <- seq_along(stage_name)

  # Build matrix for survival/maturation "season"
  sm_mat <- subdiag(n_stage, (survival * maturation)[-n_stage]) +
    diag(survival * (1 - maturation))

  # Build matrix for maternity "season"
  mat_mat <- diag(n_stage)
  mat_mat[1,] <- mat_mat[1,] + maternity

  # Combine them approriately
  if (model == "pre") {
    A <- sm_mat %*% mat_mat
  } else {
    A <- mat_mat %*% sm_mat
  }

  # Tidy up the matrix
  rownames(A) <- stage_name
  colnames(A) <- stage_name
  A <- dropzeros(A)

  A
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

#' Construct a stage-structured matrix with given mean stage durations
#'
#' @param stage_table Either a vector of stage names or a data frame with columns
#'   "stage_name", "survival", "maternity", and "duration". In the latter case the next
#'   three arguments do not need to be provided
#' @param survival A vector of stage-specific survival, on a per-timestep basis
#' @param maternity A vector of stage-specific maternities (number of offspring produced
#'   by an individual in a given stage)
#' @param duration A vector of average number of timesteps spent in each stage. If the
#'   final stage continues indefinitely (there is no maximum age) then the last element
#'   should be "Inf"
#' @param approx_method The rule for generating the matrix (see details). Defaults to
#'   "unrolled"
#' @param model Whether the matrix should represent a prebreeding ("pre") or postbreeding
#'   ("post") census model. Defaults to "post"
#' @param tol The convergence tolerance for lambda in estimating the AAS model. Defaults
#'   to 1e-6.
#'
#' @details There is no universally "best" way to construct a stage structured model based
#'   on mean stage durations. This function implements four approaches, with names based
#'   on the scheme in Kendall et al. (in review).
#'
#'   \code{approx_method = "unrolled"}: this creates an age-structured Leslie matrix where
#'   stage i is replicated \code{duration[i]} times. This is a good solution if the
#'   variance in stage duration is small, and is the best solution if there is no variance
#'   in stage duration: it correctly generates both the transient and asymptotic dynamics.
#'   The cost is a potentially large matrix, which doesn't cause R any difficulty but may
#'   be challenging to visualize (future developments will provide tools to aggregate
#'   results by stage). This method requires that all elements of \code{duration} be
#'   integers.
#'
#'   \code{approx_method = "AAS"} ("asymptotic age structure"): this creates a
#'   stage-structured Lefkovitch model where the fraction of individual maturing out of
#'   stage i is given by the formula in Crowder et al. (1994). This reproduces the mean
#'   stage durations that would be observed if the population were at the asymptotic stage
#'   structure, and (with suitable care in the sensitivity calculations) can reproduce the
#'   asymptotic calculations of the unrolled matrix. However, the model will likely not
#'   give useful estimates of transient dynamics, and will not give reliable statistics
#'   that are based on folloing cohorts (such as R0), for which SAS is preferred. As with
#'   SAS and FAS, the variance of the stage durations cannot be set, being equal to the
#'   means. This method requires an iterative approach, as the maturation fraction depends
#'   on the asymptotic growth rate, lambda1; the argument \code{tol} is used to determine
#'   when convergence has occurred. Smaller values will produce more precise estimates of
#'   the transition rates, but will require longer to calculate.
#'
#'   \code{approx_method = "SAS"} ("stable age structure"): this creates a
#'   stage-structured Lefkovitch model where the fraction of individual maturing out of
#'   stage i is given by the formula in Crouse et al. (1987). This reproduces the mean
#'   stage durations that would be observed if one is following a cohort through time, and
#'   hence is probably the best stage-based solution for calculating quantities such as
#'   R0, mean age at maturity, or some definitions of generation time. However, unless
#'   lambda = 1, this model will not give the correct mean stage duration under asymptotic
#'   conditions. As with AAS and FAS, the variance of the stage durations cannot be set,
#'   being equal to the means.
#'
#'   \code{approx_method = "FAS"} ("flat age structure"): this creates a stage-structured
#'   Lefkovitch model where the fraction of individual maturing out of stage i is
#'   \code{1/duration[i]}. This is not likely to ever be a good solution, as the mean
#'   stage duration will only be matched under very special conditions and the stage
#'   duration variance is uncontrolled (it is equal to the mean). However, it is the first
#'   solution proposed by Caswell (2001) and, being easy to calculate by hand, is popular
#'   in the literature.
#'
#' @return A projection matrix. For the "unrolled" case the row and column names are a
#'   concatenation of the stage name and the age; for the others they are the stage names.
#'
#' @references Caswell, H. 2001. Matrix population models: Construction, analysis, and
#'   interpretation, 2nd edition. Sinauer Associates, Sunderland, MA.
#'
#'   Crouse, D. T., L. B. Crowder, and H. Caswell. 1987. A stage-based population model
#'   for loggerhead sea turtles and implications for conservation. Ecology 68:1412-1423.
#'
#'   Crowder, L. B., D. T. Crouse, S. S. Heppell, and T. H. Martin. 1994. Predicting the
#'   impact of turtle excluder devices on loggerhead sea turtle populations. Ecological
#'   Applications 4:437-445.
#'
#'   Kendall, B.E., M. Fujiwara, J. Diaz-Lopez, S. Schneider, J. Voigt, and S. Wiesner. In
#'   review. Persistent problems in the construction of matrix population models.
#'   Ecological Modelling.
#'
#' @export
#'
#' @examples
make_stage4age_matrix <- function(stage_table, survival = stage_table$survival,
  maternity = stage_table$maternity, duration = stage_table$duration,
  approx_method = c("unrolled", "AAS", "SAS", "FAS"), model = c("post", "pre"),
  tol = 10^(-6)) {
  # Argument unpacking
  if(is.data.frame(stage_table)) {
    stage_name <- stage_table$stage
  } else {
    stage_name <- stage_table
  }
  approx_method <- match.arg(approx_method)
  model <- match.arg(model)

  # Useful indicies
  n_stage <- length(stage_name)
  stage_index <- seq_along(stage_name)

  if (approx_method == "unrolled") {
    if (is.infinite(duration[n_stage])) duration[n_stage] <- 1
    stage2age <- rep(stage_index, duration)
    ages <- seq_along(stage2age) - 1
    lifetable <- data.frame(x = ages,
                            sx = survival[stage2age],
                            mx = maternity[stage2age]
                            )
    A <- make_Leslie_matrix(lifetable, model = model)

    # Construct row and column names
    classnames <- paste(stage_name[stage2age], ages, sep = "")
    if (model == "pre") classnames <- classnames[-1]
    rownames(A) <- classnames
    colnames(A) <- classnames
  } else if (approx_method == "AAS") {
    lambda_old <- 0
    lambda_new <- 1
    while(abs(lambda_new - lambda_old) > tol) {
      lambda_old <- lambda_new
      surv_lam <- survival/lambda_old
      maturation <- surv_lam^(duration - 1) * (1 - surv_lam) / (1 - surv_lam^duration)
      if (is.infinite(duration[n_stage])) maturation[n_stage] <- 0
      st_table <- data.frame(stage = stage_name,
                             survival = survival,
                             maternity = maternity,
                             maturation = maturation)
      A_try <- make_Lefkovitch_matrix(st_table, model = model)
      lambda_new <- lambda1(A_try)
    }
    A <- A_try
  } else if (approx_method == "SAS") {
    maturation <- survival^(duration - 1) * (1 - survival) / (1 - survival^duration)
    if (is.infinite(duration[n_stage])) maturation[n_stage] <- 0
    st_table <- data.frame(stage = stage_name,
                           survival = survival,
                           maternity = maternity,
                           maturation = maturation)
    A <- make_Lefkovitch_matrix(st_table, model = model)
  } else if (approx_method == "FAS") {
    maturation <- 1/duration
    if (is.infinite(duration[n_stage])) maturation[n_stage] <- 0
    st_table <- data.frame(stage = stage_name,
                           survival = survival,
                           maternity = maternity,
                           maturation = maturation)
    A <- make_Lefkovitch_matrix(st_table, model = model)
  }

  A

}