## Test environments
* local OS X install, R 3.4.3
* ubuntu 12.04 (on travis-ci), R 3.4.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Sebastian J. Teran Hidalgo <sebastianteranhidalgo@gmail.com>'

This is a new submission.

* checking DESCRIPTION meta-information ... WARNING 
  Dependence on R version '3.4.2' not with patchlevel 0
    
I changed R (>= 3.4.2) to R (>= 3.4) in the DESCRIPTION.

* checking top-level files ... WARNING
Conversion of 'README.md' failed:
pandoc.exe: Could not fetch ancut.png
  
I have simplified the README file and now it runs.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

* I have run R CMD check on the NUMBER downstream dependencies.
  
* FAILURE SUMMARY

* All revdep maintainers were notified of the release on RELEASE DATE.
