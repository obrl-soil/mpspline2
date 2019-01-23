
<!-- README.md is generated from README.Rmd. Please edit that file -->
mpspline2
=========

This package is a space to re-work `GSIF::mpspline()`, because the doco says it needs tidying up.

This is only set up as a package so I could use `testthat` and `roxygen2` easily. Also I don't understand SCM repositories, and don't really want to &gt;.&gt;

Installation
------------

Don't! Or use `devtools::install_github('obrl-soil/mpspline2)`, I'm not the boss of you.

Outcome
-------

Better:

-   Full unit testing coverage, for great confidence
-   Clearer documentation and commentary (or at least, more of it :P)
-   Handles multiple input formats - don't need to use SoilProfileCollection
-   TMSE in output object

Worse:

-   Not as fast :(

???:

-   Refuses to predict 1cm values outside input data range (on purpose).
-   Default return object is different. 'Classic' option is available.

------------------------------------------------------------------------
