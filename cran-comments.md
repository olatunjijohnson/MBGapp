## Submission

This is a new submission.

## Test environments

* local: Ubuntu 22.04, R 4.5.2 (R CMD check --as-cran)
* win-builder, R-devel (Windows Server 2022): Status 1 NOTE (see below)
* mac-builder, R release (macOS): Status OK

## R CMD check results

0 errors | 0 warnings | 1 NOTE (new submission; details below).

## Notes

* This is a new release.
* "Possibly misspelled words in DESCRIPTION": these are all spelled correctly.
  "Diggle" and "Giorgi" are author surnames in the cited reference;
  "geostatistics", "geostatistical" and "Geostatistics" are standard
  terminology in the field.
* The package ships an interactive Shiny application in `inst/MBGapp`. The
  packages listed in `Imports` are required by that application and are loaded
  by it at runtime via `library()`. They are intentionally not imported into
  the package namespace (importing them would load optional system dependencies
  such as Tcl/Tk via `geoR` at package-load time). A NOTE of the form
  "Namespaces in Imports field not imported from" may therefore appear on some
  platforms. They are declared in `Imports` (rather than `Suggests`) so that
  installing the package installs everything the bundled application needs.
* `INLA` (in `Suggests`, available from the `Additional_repositories`) is an
  optional fast-inference backend used only by the bundled application. The
  application checks for it at runtime and degrades gracefully when it is not
  installed.
