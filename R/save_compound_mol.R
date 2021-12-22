#' Download compound .mol file
#'
#' Using KEGG ID values for compound and pathway, download
#' and clean up corresponding mol file
#'
#' @param compound_id character string that is 5 digits prepended with a "C"
#' @param pathway_id character string that is 5 digits prepended with "map"
#' @param path relative path to location to download data
#' @param redownload download file again even if it has already been downloaded at path
#'
#' @return downloaded .mol file
#' @export
#' @importFrom magrittr %>%
#' @importFrom utils write.table
save_compound_mol <- function(compound_id, pathway_id, path, redownload = FALSE) {
  . <- NULL
  if (!stringr::str_detect(compound_id, "^[C][:digit:]{5}$")) {
    stop("compound_id is not in the correct KEGG format")
  }
  if (!stringr::str_detect(pathway_id, "^[m][a][p][:digit:]{5}$")) {
    stop("pathway_id is not in the correct KEGG format")
  }
  pathway_dir <- paste0(path, "/", pathway_id)
  if (!dir.exists(pathway_dir)) {
    dir.create(pathway_dir, recursive = TRUE)
  }
  if (!file.exists(paste0(pathway_dir, "/", compound_id, ".mol")) | isTRUE(redownload)) {
    mol <- KEGGREST::keggGet(compound_id, option = "mol")
    mol_clean <- gsub(">.*", "", mol) %>%
      paste0("\n\n\n", .)
    file_path <- paste0(pathway_dir, "/", compound_id, ".mol")
    write.table(mol_clean,
      file = file_path, row.names = FALSE,
      col.names = FALSE, quote = FALSE
    )
  }
}

stringr::str_detect("C12345687", "[C][:digit:]{5}")
stringr::str_detect("C12345fjdklsfjdsl", "[C][0-9]{5}")
