# run_app.R
#
#' @title A function to run the app
#' @description Run the app
#' @details The function outputs the app for the user to supply the data and analysis the data
#' Runs the MBGapp Shiny web application.
#' @export
run_app <- function() {
  appDir <- system.file('MBGapp', package='MBGapp')
  if (appDir == "") {
    stop("Could not find MBGapp. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
