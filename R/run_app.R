#' Runs the SpatialEpiApp Shiny web application.
#' @export
run_app <- function() {
  appDir <- system.file('MBGapp', package='MBGapp')
  if (appDir == "") {
    stop("Could not find MBGapp. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
