## code to prepare the `teasel` dataset
## All data comes from Table 3 of Werner & Caswell (1977)
## Using "Field A"

stages <- c("seeds", "dormant1", "dormant2",
            "sm_rosette", "med_rosette", "lg_rosette", "flowering")

# Reproduction
seed_production <- c(0, 0, 0, 0, 0, 0, 431)

# Survival, based on summing the published transitions
surv <- c(0.748 + 0.008 + 0.07 + 0.002,
          0.966 + 0.013 + 0.007 + 0.008,
          0.01,
          0.125 + 0.125 + 0.038,
          0.238 + 0.245 + 0.023,
          0.167 + 0.75,
          0)

# Maturation/growth matrix, normalized so columns sum to one
Tmat <- matrix(0, 7, 7)
Tmat[2:6, 1] <- c(0.748, 0, 0.008, 0.07, 0.002) / surv[1]
Tmat[3:6, 2] <- c(0.966, 0.013, 0.007, 0.008) / surv[2]
Tmat[4, 3] <- 1
Tmat[4:7, 4] <- c(0.125, 0.125, 0.038, 0) / surv[4]
Tmat[5:7, 5] <- c(0.238, 0.245, 0.023) / surv[5]
Tmat[6:7, 6] <- c(0.167, 0.75) / surv[6]

# Create the data
teasel <- list(
  demog <- data.frame(
    stage = stages,
    seed_production = seed_production,
    survival = surv
  ),
  transitions = Tmat
)

usethis::use_data(teasel, overwrite = TRUE)
