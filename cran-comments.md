## Submission

This is a maintenance update (0.2.6 -> 0.2.7). It fixes several bugs in the
sequential test p-value adjustment (`pSeqStop()`, `gpdSeqTests()`,
`gevrSeqTests()`) and in model fitting (`gpdFit()`, `gevrFit()`), adds formal
package citations, and makes the vignette robust to `SpatialExtremes` being
unavailable. See NEWS.md for the full list of changes.

## Test environments

* local Windows 11, R 4.6.0
* win-builder (R-devel, 2026-06-19 r90183)

## R CMD check results

0 errors | 0 warnings | 0 notes

`SpatialExtremes` is used only in the vignette, is listed in `Suggests`, and all
uses are guarded by `requireNamespace("SpatialExtremes", quietly = TRUE)`, so the
package builds and checks normally whether or not it is installed.

## Reverse dependencies

We checked 3 reverse dependencies (2 from CRAN + 1 from Bioconductor), comparing
R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
