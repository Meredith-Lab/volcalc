#' Example .mol files
#' 
#' `volcalc` comes bundled with some example .mol files in its `inst/extdata`
#' directory.  This function provides easy access to them.
#'
#' @details File names are the KEGG identifiers.  Compound names are as follows:
#' - C00031: D-Glucose
#' - C00157: Phosphatidylcholine
#' - C08491: (-)-Jasmonic acid
#' - C16181: beta-2,3,4,5,6-Pentachlorocyclohexanol
#' - C16286: Geosmin
#' - C16521: Isoprene

#' @returns File paths to installed example .mol files.
#' @export
#'
#' @examples
#' #return paths to all example .mol files
#' mol_example()
#' 
#' #examine the contents of a file
#' readLines(mol_example()[1])
mol_example <- function() {
  dir(system.file("extdata", package = "volcalc"), full.names = TRUE)
}

