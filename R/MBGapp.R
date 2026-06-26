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
#' @keywords internal
#'
#' @importFrom shiny runApp
#' @importFrom geoR variog
#' @importFrom ggplot2 ggplot
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @importFrom readr read_csv
#' @importFrom tidyr pivot_longer
#' @importFrom sf st_as_sf
#' @importFrom leaflet leaflet
#' @importFrom leafem addStarsImage
#' @importFrom tidyterra geom_spatraster
#' @importFrom stars st_as_stars
#' @importFrom RiskMap glgpm
#' @importFrom terra rast
#' @importFrom grDevices colorRampPalette
#' @importFrom shinyjs useShinyjs
#' @importFrom splines ns
#' @importFrom httr2 request
#' @importFrom rmarkdown render
"_PACKAGE"
