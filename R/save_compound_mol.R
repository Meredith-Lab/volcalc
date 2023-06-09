#' Download compound .mol file
#'
#' Using KEGG ID values for compound and pathway, download and clean up
#' corresponding mol file
#'
#' @param compound_id A character string that is 5 digits prepended with a "C".
#' @param compound_formula A character string detailing a compound formula.
#' @param pathway_id An optional character string specifying KEGG pathway ID, in
#'   format of 5 digits prepended with "map".
#' @param path An optional parameter to set relative path to location to
#'   download data.
#' @param redownload Download file again even if it has already been downloaded
#'   at path.
#'
#' @return Downloaded .mol file for compound in path folder.
#'
#' @examples
#' \dontrun{
#' save_compound_mol(compound_id = "C16181")
#' }
#'
#' @export
save_compound_mol <-
  function(compound_id = NULL,
           compound_formula = NULL,
           pathway_id = NULL,
           path = "data",
           redownload = FALSE) {
  
  #assign variables to quiet devtools::check()
  . <- kegg_id <- NULL
  
  if (is.null(compound_id) & is.null(compound_formula)) {
    stop("either compound_id or compound_formula needs to be specified")
  }
  if (!is.null(compound_formula)) {
    formula_options <- KEGGREST::keggFind("compound", compound_formula, "formula") %>%
      data.frame() %>%
      tibble::rownames_to_column("kegg_id") %>%
      dplyr::mutate(kegg_id = substr(kegg_id, 5, nchar(kegg_id))) %>%
      dplyr::select(kegg_id)
    if (nrow(formula_options) == 0) {
      stop("compound_formula does not exist")
    } else if (nrow(formula_options) == 1) {
      compound_id <- formula_options$kegg_id[1]
      message("single compound found for formula, will use KEGG ID ", compound_id)
    } else {
      stop("multiple compounds found for formula, choose from these options and rerun: ", paste(as.character(unique(formula_options$kegg_id)), collapse = ", "))
    }
  }
  if (!stringr::str_detect(compound_id, "^[C][:digit:]{5}$")) {
    stop("compound_id is not in the correct KEGG format")
  }
  if (!is.null(pathway_id)) {
    if (!stringr::str_detect(pathway_id, "^[m][a][p][:digit:]{5}$")) {
      stop("pathway_id is not in the correct KEGG format")
    }
  }
  if(!is.null(pathway_id)){
    pathway_dir <- file.path(path, pathway_id)
  } else {
    pathway_dir <- path
  }
  if (!dir.exists(pathway_dir)) {
    dir.create(pathway_dir, recursive = TRUE)
  }
  if (!file.exists(file.path(pathway_dir, paste0(compound_id, ".mol"))) | isTRUE(redownload)) {
    mol <- KEGGREST::keggGet(compound_id, option = "mol")
    mol_clean <- gsub(">.*", "", mol) %>%
      paste0("\n\n\n", .)
    file_path <- file.path(pathway_dir, paste0(compound_id, ".mol"))
    utils::write.table(mol_clean,
      file = file_path, row.names = FALSE,
      col.names = FALSE, quote = FALSE
    )
  }
}
