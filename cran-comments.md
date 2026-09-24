## Submission

This is a minor update (1.2.0). It improves the peaks-over-threshold selection
and GPD fitting in the extreme value analysis pipeline, exposes the `shape_bnd`
argument in `TsEvaNs()`, `tsEvaSampleData()` and `tsGetPOT()`, reuses the GPD
fit computed during threshold selection in `tsEVstatistics()`, and adds
robustness in the trend threshold estimation. See NEWS.md for the full list of
changes.

## Test environments

* local: Windows, R release
* GitHub Actions: macOS (release), Windows (release),
  Ubuntu (devel, release, oldrel-1)
* win-builder: R release and R devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.
