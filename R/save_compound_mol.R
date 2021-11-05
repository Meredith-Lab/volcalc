#' Download compound .mol file
#'
#' Using KEGG ID values for compound and pathway, download
#' and clean up corresponding mol file
#'
#' @param compound_id character
#' @param pathway_id character
#'
#' @importFrom magrittr %>%
#' @importFrom utils write.table
save_compound_mol <- function(compound_id, pathway_id){
  . <- NULL
  pathway_dir <- paste0("data/compound_mols/", pathway_id)
  if(!dir.exists(pathway_dir)){
    dir.create(pathway_dir, recursive = TRUE)
  }
  mol <- KEGGREST::keggGet(compound_id, option = "mol")
  mol_clean <- gsub(">.*", "", mol) %>%
    paste0("\n\n\n", .)
  file_path <- paste0(pathway_dir, "/", compound_id, ".mol")
  write.table(mol_clean, file = file_path, row.names = FALSE,
              col.names = FALSE, quote = FALSE)
}
