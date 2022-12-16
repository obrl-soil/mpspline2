
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![R-CMD-check](https://github.com/obrl-soil/mpspline2/workflows/R-CMD-check/badge.svg)](https://github.com/obrl-soil/mpspline2/actions)
[![Coverage
status](https://codecov.io/gh/obrl-soil/mpspline2/branch/master/graph/badge.svg)](https://codecov.io/github/obrl-soil/mpspline2?branch=master)
[![CRAN](https://www.r-pkg.org/badges/version/mpspline2)](https://cran.r-project.org/package=mpspline2)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/mpspline2)](https://www.r-pkg.org/pkg/mpspline2)

# mpspline2

This package is a standalone re-implementation of
[`GSIF::mpspline()`](https://r-forge.r-project.org/scm/viewvc.php/pkg/R/mpspline.R?view=markup&revision=240&root=gsif),
which applies a mass-preserving spline to soil attributes. Splining soil
data is a safe way to make continuous down-profile estimates of
attributes measured over discrete, often discontinuous depth intervals.

## Installation

Install from CRAN with

``` r
install.packages('mpspline2')
```

Install from github with

``` r
devtools::install_github("obrl-soil/mpspline2")
```

## Example

``` r
library(mpspline2)
dat <- data.frame("SID" = c( 1,  1,  1,  1),
                   "UD" = c( 0, 20, 40, 60),
                   "LD" = c(10, 30, 50, 70),
                  "VAL" = c( 6,  4,  3, 10),
                   stringsAsFactors = FALSE)
dat
#>   SID UD LD VAL
#> 1   1  0 10   6
#> 2   1 20 30   4
#> 3   1 40 50   3
#> 4   1 60 70  10
spl_dat <- mpspline_tidy(obj = dat, var_name = 'VAL')
spl_dat$est_dcm
#>     SID UD  LD SPLINED_VALUE
#> 1.1   1  0   5      6.105160
#> 1.2   1  5  15      5.618341
#> 1.3   1 15  30      4.316876
#> 1.4   1 30  60      4.192058
#> 1.5   1 60 100      9.732400
```

### Asking for help

If you get stuck using this package, please post a question on [Stack
Overflow](https://stackoverflow.com/). This means that others can
benefit from the discussion, and more people are available to help you.
You’re welcome to ping me in a comment or on Mastodon ([@obrl_soil@mastodon.social](https://mastodon.social/@obrl_soil)) to get
my attention.

------------------------------------------------------------------------
