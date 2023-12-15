utils::globalVariables(".data")

# This only exists to silence a R CMD check note about ChemmineOB being in
# Imports but no functions being used in volcalc directly.
zzz <- function() {
  ChemmineOB::prop_OB
}


#' Extracts first element of results of atomcount() as a tibble
#' 
#' ChemmineR::atomcount() returns a list of table objects, sometimes named,
#' sometimes unnamed. Currently not taking advantage of the vectorization in
#' ChemmineR and only allowing length 1 lists here.
#'
#' @param atomcount output of ChemmineR::atomcount()
#'
#' @returns a tibble
#' @noRd
#' 
atomcount2tibble <- function(atomcount) {
  if (length(atomcount) != 1) {
    stop("Only length 1 lists are allowed")
  }
  atomcount[[1]] %>% 
    tibble::as_tibble_row() %>% 
    dplyr::mutate(dplyr::across(dplyr::everything(), as.integer))
}