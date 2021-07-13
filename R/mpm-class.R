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