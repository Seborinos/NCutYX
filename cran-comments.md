## Test environments
* local OS X install, R 3.4.3
* ubuntu 12.04 (on travis-ci), R 3.4.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Sebastian J. Teran Hidalgo <sebastianteranhidalgo@gmail.com>'

This is a new submission.

## Comments by CRAN

* Thanks, please add a reference for the method in the 'Description' field of your DESCRIPTION file in    the form
  authors (year) <doi:...>
  authors (year) <arXiv:...>
  authors (year, ISBN:...)
  with no space after 'doi:', 'arXiv:' and angle brackets for auto-linking.

I have added one of our published references in the DESCRIPTION as requested by the reviewer. We have other manuscripts related to this package that have been submitted but not yet published so they do not have dois. We talked among authors and decided not to submit anything to the arXiv. I hope it is acceptable to the reviewer only having one reference here. We will add more references in later versions of the package once they get published.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

* I have run R CMD check on the NUMBER downstream dependencies.
  
* FAILURE SUMMARY

* All revdep maintainers were notified of the release on RELEASE DATE.
