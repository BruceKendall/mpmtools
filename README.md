
<!-- README.md is generated from README.Rmd. Please edit that file -->
mpmtools: Tools for Building and Analyzing Matrix Population Models
===================================================================

**mpmtools** provides A collection of tools for building and analyzing matrix population models (MPMs). An MPM represents the demography of an age-, stage-, or size-structured population in one or more projection matrices. This package contains tools that help ease the proces of constructing MPMs of varying degrees of complexity, and provides functions for advanced model analysis. mpmtools is intended to help non-modelers more easily construct MPMs that reliably represent their data, and to perform analyses that would otherwise require programming skills.

Installation
------------

You can install mpmtools from github with:

``` r
# install.packages("devtools")
devtools::install_github("BruceKendall/mpmtools")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(mpmtools)

# Create a demography schedule, with juvenile and senescent age classes
demog_sched <- data.frame(x = 0:7,
                          sx = c(0.05, 0.2, 0.35, 0.8, 0.9, 0.9, 0.75, 0.4),
                          mx = c(0, 0, 0, 0.5, 1, 3, 3, 1.5))

# Construct a Leslie matrix from this demography schedule
A1 <- make_Leslie_matrix(demog_sched)
A1
#>      [,1] [,2]  [,3] [,4] [,5] [,6]  [,7]
#> [1,]  0.0 0.00 0.025 0.05 0.15 0.15 0.075
#> [2,]  0.2 0.00 0.000 0.00 0.00 0.00 0.000
#> [3,]  0.0 0.35 0.000 0.00 0.00 0.00 0.000
#> [4,]  0.0 0.00 0.800 0.00 0.00 0.00 0.000
#> [5,]  0.0 0.00 0.000 0.90 0.00 0.00 0.000
#> [6,]  0.0 0.00 0.000 0.00 0.90 0.00 0.000
#> [7,]  0.0 0.00 0.000 0.00 0.00 0.75 0.400

# Calculated the asymptotic growth rate of the population governed by this demography
#   schedule:
lambda1(A1)
#> [1] 0.5535262
```
