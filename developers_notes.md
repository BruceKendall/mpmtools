# Developer's notes

In general, this package is written according to the style recommended in Hadley Wickhams **R Packages** book (http://r-pkgs.had.co.nz). Exceptions:

- We use the git-flow branching model (https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow). This means that the master branch has the latest non-development version, allowing "stable releases" without pushing to CRAN (once the package is closer to complete then CRAN submission will happen).
- We use the auto-generated citation (at least until a publication is written about it). This requires manually updating the "Date/Publication" field in DESCRIPTION, which is a pain!