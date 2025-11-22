## Release Summary

Feature addition (processing multiple input variables in one call), column name bugfix, checks against newer R versions

## Test environments

  * Local: Windows 11, R 4.5.1
  * Github Actions via usethis::use_github_action('check-standard')
    * Note intermittent issue with macos build - related to github action setup stages, not the package itself
  * R-devel with devtools::check_win_devel()
  * MacOS with devtools::check_mac_release()

## R CMD Check Results

  * 0 errors | 0 warnings | 0 notes
  
## Downstream dependencies

There are currently no downstream dependencies for this package. Suggested 'aqp' will be unaffected.
