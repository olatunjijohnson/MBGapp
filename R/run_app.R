#' Launch the MBGapp Shiny application
#'
#' Runs the interactive Shiny application for model-based geostatistical
#' analysis that is bundled with the package. The application lets the user
#' upload data, explore it, fit geostatistical models and produce predictions
#' and reports.
#'
#' @return No return value; called for the side effect of launching the Shiny
#'   application.
#' @examples
#' \dontrun{
#' run_app()
#' }
#' @export
run_app <- function() {
  appDir <- system.file("MBGapp", package = "MBGapp")
  if (appDir == "") {
    stop("Could not find the MBGapp application directory. Try re-installing 'MBGapp'.",
         call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
