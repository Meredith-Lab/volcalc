#' Get path to example .mol file
#' 
#' `volcalc` comes bundled with some example .mol files in its `inst/extdata`
#' directory.  This function provides easy access to them.
#'
#' @param path Name of file.  If `NULL` (default), the options for example .mol
#'   files will be listed instead.
#'
#' @export
#'
#' @examples
#' mol_example()
#' mol_example("C16181.mol")
mol_example <- function(path = NULL) {
  if (is.null(path)) {
    dir(system.file("extdata", package = "volcalc"))
  } else {
    system.file("extdata", path, package = "volcalc", mustWork = TRUE)
  }
}

