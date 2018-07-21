context("Prebreeding to postbreeding conversion")

# Create some simple matrices to work with
# 1. The age-structured model from Kendall et al. (2019)
S0 <- 0.2
Sx <- c(0.4, 0.4, 0.9)
bx <- c(0, 0, 3)
A1_pre <- matrix(c(S0 * bx,
                   Sx[1], 0, 0,
                   0, Sx[2], Sx[3]), 3, 3, byrow = TRUE)

# 2. The stage-structured model from Kendall et al. (2019)
S0 <- 0.2
Si <- c(0.4, 0.9)
bi <- c(0, 3)
gamma_1 <- 0.5
A2_pre <- matrix(c((1 - gamma_1) * Si[1], S0 * bi[2],
                   gamma_1 * Si[1], Si[2]), 2, 2, byrow = TRUE)

test_that("pre and post have same lambda1", {
  lambda1 <- function(A) eigen(A)$values[1]
  expect_equal(lambda1(A1_pre), lambda1(suppressWarnings(pre_to_post(S0, A1_pre))))
  expect_equal(lambda1(A2_pre), lambda1(suppressWarnings(pre_to_post(S0, A2_pre))))
})

test_that("null Fmat generates warning", {
  expect_warning(pre_to_post(S0, A1_pre))
})

test_that("unsupported Fmat generates error", {
  expect_error(pre_to_post(S0, A1_pre, matrix(1, 3, 3)))
})

