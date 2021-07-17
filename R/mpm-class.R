new_mpm <- function(demog = list(births = numeric(),
                                 survival = numeric(),
                                 transitions = numeric()),
                    timestep = character(),
                    matrix_type = character(),
                    class_type = character(),
                    class_names = character(),
                    newborn_classes = character(),
                    census = character()) {
  # HOLD these until I know whether any are needed
  # stopifnot(c("births", "survival", "transitions") %in% names(demog))
  # num_class <- length(demog$survival)
  # if (is.na(class_names)) {
  #   class_names <- paste(class_type, 0:(num_class - 1))
  # }
  # stopifnot(newborn_classes %in% class_names)

  num_class <- length(demog$survival)
  nb_idx <- match(newborn_classes, class_names)

  if (matrix_type == "Leslie") {
    # Assumptions:
    #   - `births` and `survival` are vectors with length equal to the number of
    #     age classes
    #   - The first age class is newborn
    #   - If the final element of `survival` is nonzero then there is a self-
    #     loop in the final age class
    Bmat <- diag(num_class)
    Bmat[1, ] <- demog$births
    Pmat <- diag(demog$survival)
    Tmat <- subdiag(num_class, 1)
    if (demog$survival[num_class] != 0) Tmat[num_class, num_class] <- 1
  } else if (matrix_type == "Lefkovitch") {
    # Assumptions:
    #   - `births` `survival`, and 'maturation`` are vectors with length equal to
    #     the number of age classes
    #   - The first age class is newborn
    #   - maturation[i] is the probability of moving from class i to class i+1;
    #     the remainder stay in the same class
    Bmat <- diag(num_class)
    Bmat[1, ] <- demog$births
    Pmat <- diag(demog$survival)
    Tmat <- diag(1 - demog$maturation)
    Tmat <- subdiag(Tmat, demog$maturation[-num_class])
  } else if (matrix_type == "general") {
    # Assumptions:
    #   - `births` `survival`, and 'maturation`` are vectors with length equal to
    #     the number of age classes
    #   - The first age class is newborn
    #   - maturation[i] is the probability of moving from class i to class i+1;
    #     the remainder stay in the same class
    Bmat <- diag(num_class)
    if (is.vector(demog$births)) {
      Bmat[nb_idx, ] <- demog$births
    } else {
      Bmat[nb_idx, ] <- t(demog$births)
    }
    Pmat <- diag(demog$survival)
    Tmat <- demog$transitions
  } else {
    stop("Matrix type ", matrix_type, "not implemented!")
  }

  if (length(class_names == nrow(Bmat))) {
    colnames(Bmat) <- rownames(Bmat) <- class_names
    colnames(Pmat) <- rownames(Pmat) <- class_names
    colnames(Tmat) <- rownames(Tmat) <- class_names
  }
  TP <- Tmat %*% Pmat
  Apre <- dropzeros(TP %*% Bmat)
  Apost <- dropzeros(Bmat %*% TP)

  structure(list(demog = demog,
                 Apre = Apre,
                 Apost = Apost,
                 Bmat = Bmat,
                 Pmat = Pmat,
                 Tmat = Tmat,
                 class_names = class_names,
                 newborn_classes = newborn_classes),
            timestep = timestep,
            matrix_type = matrix_type,
            class_type = class_type,
            census = census,
            class = "mpm")
}

#' Construct an mpm (matrix population model) object
#'
#' This function takes information about population structure, demography (e.g.,
#' birth and survival rates), and growth/maturation, and correctly does the
#' necessary accounting to create the projection matrix (or matrices) that
#' represent the biological information. For time-invariant models, the
#' projection matrix can be extracted and analyzed directly, although the
#' package has analysis tools that take advantage of the underlying demographic
#' information that is encoded in an `mpm` object. When more complex models are
#' implemented (e.g., with time-varying parameters or density-dependence), then
#' the `mpm` object will encode the information needed to simulate or analyze
#' the model.
#'
#' The contents of `demog_infa` vary depending on which model type is requested.
#' It should be a list or data frame, with each element or column named. In all
#' models, a column or list element giving unique identifiers for each class is
#' required *unless* both `class_type` and `class_names` are specified. The list
#' element or column containing the vector of class identifiers should either
#' have the same name as is given for `class_type` (e.g., "age") or be the first
#' element of `demog_info`.
#'
#' Most (all?) models requirw that `denog_info` contain elements that specify
#' class-specific birth and survival rates. Below these are described as having
#' names "births" and "survival"; you may also use "mx" and "Px", respectively,
#' for consistency with life table notation.
#'
#' ## Model types
#'
#' \describe{
#'
#' \item{`matrix_type = "Leslie"`}{Creates a classic single-sex, age-structured
#' (Leslie) model. The only required elements of `demog_info` are `births` (a
#' vector of age-specific birth rates, including zero for newborns; equivalent
#' to \eqn{m_x} in a life table) and `survival` (a vector of age specific
#' survival rates, equivalent to \eqn{P_x} in a life table; if the last element
#' is non-zero, then the last age-class will include a self-loop).}
#'
#' \item{`matrix_type = "Lefkovitch"`}{Creates a classic single-sex,
#' stage-structured (Lefkovitch) model, in which individuals in stage x, if the
#' survive, either stay in their stage (with probability 1 - g_x) or mature
#' to the next stage (with probability g_x). The only required elements
#' of `demog_info` are `births` (a vector of stage-specific birth rates,
#' including zero for newborns); `survival` (a vector of stage specific survival
#' rates); and `maturation` (a vector of stage-specific maturation rates; note that
#' the last stage should have a maturation rate of zero).}
#'
#' }
#'
#' ## Census
#'
#' The value of `census` gives the preferred census time for representing the
#' model itself as well as certain outcomes (e.g., which stages to include in
#' model projections). For birth-pulse model, the standard options are "post" (a
#' postbreeding-census model, which looks at the population immediately after
#' births have occurred, so that the youngest individuals are newborns) and
#' "pre" (a prebreeding-census model, which looks at the population immediately
#' before births occur, so that the youngest individuals are one timestep old).
#' The choice of a census time is important for visualization, but does not
#' affect the model's dynamic properties. In addition, the value of census that
#' is encoded in the `mpm` object can always be overridden in any function in
#' which it affects the output (e.g., projection of population vectors).
#'
#' @param demog_info A list or data frame (or equivalent, e.g., a tibble)
#'   containing information on the population structure and demography. See the
#'   **Details** section for more information.
#' @param matrix_type A character string specifying the model type. Must match
#'   one of the values described in **Details**
#' @param census The moment within the timestep at which the model "censuses"
#'   the population. For non-seasonal birth-pulse models, allowable values are
#'   "pre" (for a prebreeding-census model) or "post" (for a postbreeding-census
#'   model). The choice affects structure of the matrix and the length of the
#'   population vector, but does not affect the model's dynamics or asymptotic
#'   properties (see **Details**). Defaults to "post", which includes newborns
#'   in the population vector.
#' @param timestep The timestep of the model; defaults to "year". Not currently
#'   used for anything (other than being reported by `print`), but may become
#'   useful with future implementations of seasonal and plotting methods.
#' @param class_type A conceptual descriptor of the classes that structure the
#'   population (e.g., "age", "stage", or "size"). If not specified, then it is
#'   taken from the name of the first element of `demog_info`.
#' @param class_names A character vector giving the names of the classes. If not
#'   specified, then it is constructed by concatenating `class_type` with the
#'   element of `demog_info` that lists class names of indices.
#' @param newborn_classes A character vector (often of length one) giving the
#'   names of the classes that should be treated as newborn individuals. All
#'   elements of this vector should also be in `class_names`. Defaults to the
#'   first element of `class_names`.
#'
#' @return Ao object of class `mpm`, which is a list with the following slots:
#'
#'   - `demog`: a list containing the demographic parameters extracted from `x`
#'
#'   - `Apre` and `Apost`: the prebreeding- and postbreeding-census matrices.
#'   These can be accessed easily using `Amat()`.
#'
#'   - `Bmat`, `Pmat`, and `Gmat`: the birth, survival and transition ("growth")
#'   matrices. Users do not normally need to access these; they are stored for
#'   use by the analysis functions. - `class_names` and `newborn_classes`:
#'   character vectors giving the names of all the classes and the newborn
#'   classes, respectively.
#'
#'   It also has attributes `timestep`, `matrix_type`, `class_type`, and
#'   `census`, which store the values provided in the corresponding input
#'   parameters.
#'
#'   The `print` method gives a formatted report on the metadata as well as the
#'   projection matrix associated with the value of `census`.
#' @export
#'
#' @examples
#' # Leslie matrix model
#' bx <- c(0, 0.040, 1.470, 2.065, 2.440, 3.250, 3.250, 3.250) # Births
#' px <- c(0.424, 0.726, 0.513, 0.361, 0.175, 0.700, 0.286, 0) # Survival
#' my_lifetable <- data.frame(age = 0:7, mx = bx, Px = px)
#' mpm(my_lifetable, "Leslie") # default postbreeding-census model
#' mpm(my_lifetable, "Leslie", census = "pre") # prebreeding-census model
#' # Specify class names
#' cnames <- c("newborn", as.character(1:7))
#' mpm(my_lifetable, "Leslie", class_names = cnames)
#'
#' # Make a Lefkovitch stage structured model
#' bx <- c(0, 0, 303)
#' px <- c(0.92 * 0.03 * 0.55, 0.36, 0.69)
#' gx <- c(1, 0.09/0.36, 0)
#' my_lt <- data.frame(births = bx, survival = px, maturation = gx)
#' mpm(my_lt, "Lefkovitch", class_type = "stage",
#'     class_names = c("eggs", "juveniles", "adults"))
mpm <- function(demog_info, matrix_type,
                census = "postbreeding",
                timestep = "year",
                class_type = character(),
                class_names = character(),
                newborn_classes = character()) {
  # Argument checking
  # - demog_info should be either a data frame (or tibble) or a list
  # - matrix_type and census need to match a list of possible values
  # - Lengths/dimensions of various objects should match
  if (!is.list(demog_info)) {
    stop("demog_info should be a data frame, tibble, or list")
  }

  # Allow lifetable notation
  if ("mx" %in% names(demog_info) & !("births" %in% names(demog_info))) {
    names(demog_info)[names(demog_info) == "mx"] <- "births"
  }
  if ("Px" %in% names(demog_info) & !("survival" %in% names(demog_info))) {
    names(demog_info)[names(demog_info) == "Px"] <- "survival"
  }

  if (!all(c("births", "survival") %in% names(demog_info))) {
    stop("demog_info must contain elements births and survival")
  }
  # The following line should only apply to a standard birth-pulse model!
  census <- match.arg(census, c("prebreeding", "postbreeding"))

  # If not specified, extract sensible values for class_type, class_names,
  #   and newborn_classes from information in x and matrix_type
  if (length(class_type) == 0) {
    class_type <- names(demog_info)[1]
    cat("`class_type` not specified, setting it to", class_type, "\n")
  }
  if (length(class_names) == 0) {
    if (is.character(demog_info[[class_type]])) {
      class_names <- demog_info[[class_type]]
    } else {
      class_names <- paste(class_type, demog_info[[class_type]], sep = "_")
    }
    cat("`class_names` not specified, setting them to", class_names, "\n")
  }
  if (length(newborn_classes) == 0) {
    newborn_classes <- class_names[1]
    cat("`newborn_classes` not specified, setting it to", newborn_classes, "\n")
  }
  stopifnot(all(newborn_classes %in% class_names))

  # Type-specific requirements
  if (matrix_type == "Leslie") {
    stopifnot(exprs = {
      length(demog_info$births) == length(class_names)
      length(demog_info$survival) == length(class_names)
      length(newborn_classes) == 1
      newborn_classes == class_names[1]
    })
  } else   if (matrix_type == "Lefkovitch") {
    stopifnot(exprs = {
      length(demog_info$births) == length(class_names)
      length(demog_info$survival) == length(class_names)
      length(demog_info$maturation) == length(class_names)
      length(newborn_classes) == 1
      newborn_classes == class_names[1]
    })
  } else if (matrix_type == "general") {
    stopifnot(exprs = {
      length(demog_info$births) == length(class_names)
      length(demog_info$survival) == length(class_names)
      is.matrix(demog_info$transitions)
      nrow(demog_info$transitions) == ncol(demog_info$transitions)
      all(colSums(demog_info$transitions) == 1)
    })
  }


  # Construct demog
  if (matrix_type == "Leslie") {
    demog <- list(births = demog_info$births,
                  survival = demog_info$survival)
  } else if (matrix_type == "Lefkovitch") {
    demog <- list(births = demog_info$births,
                  survival = demog_info$survival,
                  maturation = demog_info$maturation)
  } else if (matrix_type == "general") {
    demog <- list(births = demog_info$births,
                  survival = demog_info$survival,
                  transitions = demog_info$transitions)
  } else {
    stop("Matrix type ", matrix_type, "not implemented!")
  }

  # Construct the mpm object
  new_mpm(demog, timestep, matrix_type, class_type, class_names,
          newborn_classes, census)
}

#' Extract a projection matrix from an `mpm` object
#'
#' `Amat` extracts the projection matrix (often called the "A matrix") from a
#' matrix population model (`mpm`) object. By default, the census coded into the
#' mpm object is used, but this can be overridden by specifying a value for
#' `census`.
#'
#' @param x An `mpm` object.
#' @param census A string specifying the census to use. Defaults to the value
#'   encoded in `x`.
#'
#' @return A projection matrix, which is a matrix object with row and column
#'   names.
#' @export
#'
#' @examples
#' # Make a Leslie matrix model with postbreeding census
#' bx <- c(0, 0.040, 1.470, 2.065, 2.440, 3.250, 3.250, 3.250)     # Births
#' px <- c(0.424, 0.726, 0.513, 0.361, 0.175, 0.700, 0.286, 0)     # Survival
#' my_lifetable <- data.frame(age = 0:7, mx = bx, Px = px)
#' my_mpm <- mpm(my_lifetable, "Leslie")
#'
#' # Extract the projection matrices
#' Amat(my_mpm)            # Postbreeding census
#' Amat(my_mpm, "pre")     # Prebreeding census
Amat <- function(x, census = attr(x, "census")) {
  if (!is.na(pmatch(census, "prebreeding"))) {
    A <- x$Apre
  } else if (!is.na(pmatch(census, "postbreeding"))) {
    A <- x$Apost
  } else {
    stop("No method implemented for census type ", census)
  }
  A
}

#' @export
print.mpm <- function(x, census = attr(x, "census"), digits = 3L, ...) {
  # Expand census
  if (!is.na(pmatch(census, c("postbreeding", "prebreeding")))) {
    census <- match.arg(census, c("postbreeding", "prebreeding"))
  }

  # Print metadata
  cat("\n", iart(attr(x, "class_type"), TRUE), "-structured ",
      attr(x, "matrix_type"), " matrix population model\n", sep = "")
  cat("    with ", iart(census), " census and ", iart(attr(x, "timestep")),
      " time step:\n\n", sep = "")

  # Make and print a pretty version of the A matrix
  A <- Amat(x, census)
  A <- matrix(as.character(round(A, digits)), nrow(A), dimnames = dimnames(A))
  print(A, quote = FALSE)
  invisible(x)
}