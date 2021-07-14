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

mpm <- function(x, matrix_type,
                census = "post",
                timestep = "year",
                class_type = character(),
                class_names = character(),
                newborn_classes = character()) {
  # Argument checking
  # - x should be either a data frame (or tibble) or a list
  # - matrix_type and census need to match a list of possible values
  # - Lengths/dimensions of various objects should match
  if (!is.list(x)) {
    stop("x should be a data frame, tibble, or list")
  }
  if (!all(c("mx", "Px") %in% names(x))) {
    stop("x must contain elements mx and Px")
  }
  # The following line should only apply to a standard birth-pulse model!
  census <- match.arg(census, c("prebreeding", "postbreeding"))

  # If not specified, extract sensible values for class_type, class_names,
  #   and newborn_classes from information in x and matrix_type
  if (length(class_type) == 0) {
    class_type <- names(x)[1]
    cat("`class_type` not specified, setting it to", class_type, "\n")
  }
  if (length(class_names) == 0) {
    if (is.character(x[[class_type]])) {
      class_names <- x[[class_type]]
    } else {
      class_names <- paste(class_type, x[[class_type]], sep = "_")
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
      length(x$mx) == length(class_names)
      length(x$Px) == length(class_names)
      length(newborn_classes) == 1
    })
  }

  # Construct demog
  if (matrix_type == "Leslie") {
    demog <- list(births = x$mx,
                  survival = x$Px)
  } else {
    stop("Matrix type ", matrix_type, "not implemented!")
  }

  # Construct the mpm object
  new_mpm(demog, timestep, matrix_type, class_type, class_names,
          newborn_classes, census)
}