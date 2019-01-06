### data.R
### Documentation for datasets

#' Example 2.1 from Caswell (2001)
#'
#' A 3 x 3 Leslie matrix model. Given as example 2.1, equation (2.8).
#'
#' @format A matrix with 3 rows and 3 columns
#'
#' @source Caswell, H. 2001. Matrix population models: Construction, analysis, and
#' interpretation. 2nd edition. Sinauer Associates, Sunderland, MA.
"Caswell_Ex_2.1"

#' Loggerhead sea turtle stage-based life tables
#'
#' Two stage-based life tables (differing in how adults are broken out) of the Little Cumberland Island loggerhead sea turtle
#' population, set up to use as
#' input to \code{\link{make_stage4age_matrix}}. Note that these tables do not exactly
#' generate the stage-based matrices published in the sources. The matrix in Table 4 of
#' Crouse et al. (1987) has inappropriate reproduction matrix elements for a
#' postbreeding-census matrix, and has an apparent typo in the survival of adults. The
#' matrix in Table 2 of Crowder et al. (1994) differs in the third decimal place of
#' subadult reproduction, probably as a consequence of rounding the maturation fraction
#' before multiplying by the adult maternity.
#'
#' @format A list with two elements, each being a data frame with four named columns.
#'
#' @source \code{loggerhead[[1]]} is derived from Table 3 of Crouse et al. (1987), and
#'   \code{loggerhead[[2]]} is Table 1 of Crowder et al. (1994).
#'
#' @references Crouse, D. T., L. B. Crowder, and H. Caswell. 1987. A stage-based
#'   population model for loggerhead sea turtles and implications for conservation.
#'   Ecology 68:1412-1423.
#'
#'   Crowder, L. B., D. T. Crouse, S. S. Heppell, and T. H. Martin. 1994. Predicting the
#'   impact of turtle excluder devices on loggerhead sea turtle populations. Ecological
#'   Applications 4:437-445.
"loggerhead"