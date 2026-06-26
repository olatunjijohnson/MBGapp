## Submission

This is a new submission.

## Test environments

* local: Ubuntu 22.04, R 4.5.2
* win-builder (devel and release) — to be run
* R-hub (Windows, macOS, Linux) — to be run

## R CMD check results

0 errors | 0 warnings | notes as below.

## Notes

* This is a new release.
* The package ships an interactive Shiny application in `inst/MBGapp`. The
  packages listed in `Imports` are required by that application at runtime and
  are imported into the package namespace so the application has access to them.
* `INLA` (in `Suggests`, available from the `Additional_repositories`) is an
  optional fast-inference backend. The application checks for it at runtime and
  degrades gracefully when it is not installed.
