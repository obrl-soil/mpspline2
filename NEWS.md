# v. 0.1.6

  * `mpspline_datchk()` bugfix per PR 2


# v. 0.1.5 [CRAN]

  * `mpspline_datchk()` now removes 0-thickness horizons (where upper and lower depths are identical)
  * `mpspline2_datchk()` now removes horizons where upper depth is greater than lower depth
  * `mpspline_datchk()` now returns messages about each edit made to the input data

# v. 0.1.4

  * bugfix for `mpspline_tidy()`, now handles input profile ID column names correctly.

# v. 0.1.3

  * output styles other than default shifted to wrapper functions
  * revised RMSE calculation from BM (resolves issue #4)
  * all component functions operate on a per-profile basis

# v. 0.1.2

  * Added option for SoilProfileCollection-style output from `mpspline()`.

# v. 0.1.1

  * Internal: correct S4 access methods for SoilProfileCollection objects (h/t Dylan Beaudette)
  * Outputs in 'default' mode now have `names()` attributes where they refer to a range e.g. "000_005cm"

# v. 0.1.0

  * Added a `NEWS.md` file to track changes to the package.
  * And did everything else.
