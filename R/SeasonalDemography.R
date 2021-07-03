#' Create a representation of structured demography over a single season
#'
#' The `SeasonalDemography` class contains information about the survival,
#' maternity, and class transitions (aging, growth, maturation) of a structured
#' population over a single season. The season may be an actual part of the
#' calendar, or it may be a conceptual part of the annual cycle (e.g.,
#' "breeding" or "nonbreeding" for a simple birth-pulse model, or "northward
#' migration" for a migratory population). It is a building block of the
#' `Demography` class, and many users will not need to access it directly.
#'
#' The constructor function requires that all slots be specified, and performs
#' basic validity checking on dimensions and vector lengths.
#'
#' At present, there are no functions for extracting or replacing the values in
#' indvidual slots in an existing `SeasonalDemography` object. The slots can, of
#' course, be accessed directly, but there will be no validity checking.
#'
#' @slot class_names A character vector (or something that can be coerced to
#'   one) giving the names of classes (ages, stages, etc.) of the population.
#'   Age classes may be provided as a numeric vector; it will be silently
#'   converted to character.
#' @slot season The name of the season. May be any text string.
#' @slot class_type The type of structuring variable (e.g., "age", "stage",
#'   "age-size"). May be any text string.
#' @slot newborn_classes The name(s) of the class(es) that represent newborn
#'   individuals. If there is only one such class, a string; otherwose a
#'   character vector. All elements of `newborn_classes` must be found in
#'   `class_names`.
#' @slot survival A numeric vector giving the survival probabilities of each
#'   class over the entire season.
#' @slot maternity A numeric matrix giving the number of newborns in each class
#'   produced by each of the non-newborn classes over the course of the season.
#'   There should be a row for each newborn class and a column for each
#'   non-newborn class. Newborns are assumed to have zero maternity; if
#'   `ncol(maternity) == length(class_names)`, then the newborn classes will be
#'   dropped. If there is a single newborn class, `maternity` may be provided as
#'   a numeric vector.
#' @slot transition A numeric matrix giving the class transitions. The value in
#'   `transition[i, j]` gives the probability (conditional on survival) that an
#'   individual starting the season in class `j` ends the season in class `i`.
#'   Because the transitions are conditional on survival, columns should sum to
#'   one.
#'
#' @return
#' An object of class `SeasonalDemogaphy`
#'
#' @export SeasonalDemography
#'
#' @examples
#' SeasonalDemography(c("Newborn", "Juvenile", "Adult"),
#'                    season = "nonbreeding",
#'                    class_type = "stage",
#'                    newborn_classes = "Newborn",
#'                    survival = c(0.1, 0.5, 0.9),
#'                    maternity = rep(0, 3),
#'                    transition = matrix(c(0, 0, 0, 1, 0.8, 0, 0, 0.2, 1),
#'                                        nrow = 3, byrow = TRUE))
setClass(
  "SeasonalDemography",
  slots = c(
    class_names = "character",
    season = "character",
    class_type = "character",
    newborn_classes = "character",
    survival = "numeric",
    maternity = "matrix",
    transition = "matrix"
  )
)

#' @rdname SeasonalDemography-class
SeasonalDemography <- function(class_names, season, class_type, newborn_classes,
                               survival, maternity, transition) {

  num_class <- length(class_names)
  num_nb <- length(newborn_classes)

  # Allow class names to be passed in as numeric vectors
  class_names <- as.character(class_names)
  newborn_classes <- as.character(newborn_classes)

  # Check for proper newborn class names
  nb_idx <- which(class_names %in% newborn_classes)
  if (length(nb_idx) != num_nb) {
    stop("One or more of the newborn class names ",
      "is not in the list of all class names")
  }

  # Get @maternity into the right shape
  if (is.null(dim(maternity))) {
    maternity <- matrix(maternity, nrow = 1)
  }
  if(ncol(maternity) == num_class) {
    maternity <- maternity[, -nb_idx, drop = FALSE]
  }

  # Confirm dimensions of @maternity, @survival and @transition before
  #   applying names
  if(ncol(maternity) != (num_class - num_nb)) {
    stop("The length or number of columns of @maternity must equal either \n",
         "the number of class names or the number of non-newborn class names.")
  }
  if(nrow(maternity) != num_nb & num_nb > 1) {
    stop("The number of rows of @maternity must equal the number of ",
         "newborn class names.")
  }
  if (length(survival) != num_class) {
    stop("@survival and @class_names must have the same length.")
  }
  if (!identical(dim(transition), c(num_class, num_class))) {
    stop("The number of rows and columns of @transition must equal the ",
         "number of class names.")
  }

  # Add names to the vectors and matrices
  names(survival) <- class_names
  rownames(transition) <- class_names
  colnames(transition) <- class_names
  colnames(maternity) <- class_names[-nb_idx]
  rownames(maternity) <- newborn_classes

  new("SeasonalDemography",
      class_names = class_names,
      season = season,
      class_type = class_type,
      newborn_classes = newborn_classes,
      survival = survival,
      maternity = maternity,
      transition = transition)
}

setValidity("SeasonalDemography", function(object) {
  if (min(object@newborn_classes %in% object@class_names) == 0) {
    "All elements of @newborn_classes must be in @class_names"
  } else if (length(object@survival) != length(object@class_names)) {
    "@survival and @class_names must be same length"
  } else if (!is.matrix(object@maternity)) {
    "@maternity must be a matrix"
  } else if (nrow(object@maternity) != length(object@newborn_classes)) {
    "@maternity must have a row for each element of @newborn_classes"
  } else if (ncol(object@maternity) != length(object@class_names) -
             length(object@newborn_classes)) {
    "@maternity must have a column for each non-newborn element of @class_names"
  } else if (!is.matrix(object@transition)) {
    "@transition must be a matrix"
  } else if (nrow(object@transition) != ncol(object@transition)) {
    "@transition must be a square matrix"
  } else if (nrow(object@transition) != length(object@class_names)) {
    "The rank of @transition must equal the length of @class_names"
  } else if (!is.numeric(object@maternity)) {
    "The elements of @maternity must be numeric"
  } else if (!is.numeric(object@transition)) {
    "The elements of @transition must be numeric"
  } else if (!all(colSums(object@transition) == 1)) {
    "The columns of @transition must each sum to one"
  } else {
    TRUE
  }
})

setMethod("show", "SeasonalDemography", function(object){
  cat(stringr::str_to_sentence(object@class_type),
      "-structured demography over the ", object@season, " season\n\n",
      "Survival rates: \n",
      sep = "")
  print(object@survival)
  cat("\nMaternity rates: \n")
  print(object@maternity)
  cat("\nTransition rates: \n")
  print(object@transition)
})