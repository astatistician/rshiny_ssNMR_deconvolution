#' @export
run_ssNMR <- function(){
	
	appDir <- system.file("app", package = "ssNMR")
	shiny::runApp(appDir)
	
}



