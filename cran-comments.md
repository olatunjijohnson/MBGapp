## Submission

This is a new submission.

## Test environments

* local: Ubuntu 22.04, R 4.5.2
* mac-builder (macos, R release)
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | notes as below.

## Notes

* This is a new release.
* The package ships an interactive Shiny application in `inst/MBGapp`. The
  packages listed in `Imports` are required by that application and are loaded
  by it at runtime via `library()`. They are intentionally not imported into
  the package namespace, so a check NOTE of the form "Namespaces in Imports
  field not imported from" is expected. They are declared in `Imports` (rather
  than `Suggests`) so that installing the package installs everything the
  bundled application needs.
* `INLA` (in `Suggests`, available from the `Additional_repositories`) is an
  optional fast-inference backend used only by the bundled application. The
  application checks for it at runtime and degrades gracefully when it is not
  installed.
