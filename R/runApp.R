runApp <- function() {
  appDir <- system.file("ProteATO", package = "ProteATO")
  if (appDir == "") {
    stop("Could not find Shiny app directory. Try re-installing the package.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}