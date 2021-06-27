# Developer's notes

In general, this package is written according to the style recommended in Hadley Wickhams **R Packages** book (http://r-pkgs.had.co.nz). Exceptions:

- We use the git-flow branching model. This means that the master branch has the latest non-development version, allowing "stable releases" without pushing to CRAN (once the package is closer to complete then CRAN submission will happen).
  - Be sure that the git-flow git library is installed (https://github.com/nvie/gitflow/wiki/Installation)
  - Useful references:
    - https://jeffkreeftmeijer.com/git-flow/
    - https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow
- We use the auto-generated citation (at least until a publication is written about it). This requires manually updating the "Date/Publication" field in DESCRIPTION, which is a pain!

## Handy reminders
- To install vignettes (to see how they look in the help system) use `install(build_vignettes = TRUE)`