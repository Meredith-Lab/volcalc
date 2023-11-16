#' Get path to example .mol file
#' 
#' `volcalc` comes bundled with some example .mol files in its `inst/extdata`
#' directory.  This function provides easy access to them.
#'
#' @param path Name of file.  If `NULL` (default), the options for example .mol
#'   files will be listed instead.
#' @details File names are the KEGG identifiers.  Compound names are as follows:
#' - C00031: D-Glucose
#' - C08491: (-)-Jasmonic acid
#' - C16181: beta-2,3,4,5,6-Pentachlorocyclohexanol
#' - C16286: Geosmin
#' - C16521: Isoprene
#'
#' @export
#'
#' @examples
#' #list all examples
#' mol_example()
#' 
#' #return path to specific example file
#' mol_example("C16181.mol")
mol_example <- function(path = NULL) {
  if (is.null(path)) {
    dir(system.file("extdata", package = "volcalc"))
  } else {
    system.file("extdata", path, package = "volcalc", mustWork = TRUE)
  }
}

