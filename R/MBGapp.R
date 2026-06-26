#' MBGapp: Interactive Shiny Application for Model-Based Geostatistics
#'
#' The MBGapp package provides a Shiny web application for teaching and applied
#' geostatistical analysis. Users can explore spatial data interactively, assess
#' spatial correlation through the empirical variogram, fit model-based
#' geostatistical models for continuous, prevalence and count outcomes, generate
#' spatial predictions, and download reports.
#'
#' The application is launched with \code{\link{run_app}}.
#'
#' The packages listed in the \code{Imports} field are required by the bundled
#' Shiny application (in \code{inst/MBGapp}) and are loaded by it at runtime;
#' they are intentionally not imported into the package namespace so that the
#' package itself loads quickly and without pulling in optional system
#' dependencies (e.g. Tcl/Tk via \pkg{geoR}) on headless machines.
#'
#' @keywords internal
"_PACKAGE"
