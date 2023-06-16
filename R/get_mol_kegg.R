#' Download compound .mol files from KEGG
#'
#' Downloads mol files corresponding to individual compounds or compounds in a
#' pathway from KEGG.
#'
#' @param compound_ids character vector of KEGG compound IDs---5 digits prepended with a "C".
#' @param pathway_ids character vector of KEGG pathway IDs---5 digits prepended with "map".
#' @param dir path to a folder to save mol files in.
#'
#' @returns a tibble with the columns `compound_ids`, `pathway_ids` (if used),
#'   and `mol_paths` (paths to downloaded .mol files)
#' @export
#'
#' @examples
#' \dontrun{
#' get_mol_kegg(compound_ids = c("C16181", "C06074"))
#' get_mol_kegg(pathway_ids = "map00361")
#' }
get_mol_kegg <- function(compound_ids, pathway_ids, dir){
  
  if(!fs::dir_exists(dir)) {
    fs::dir_create(dir)
  }
  #TODO check inputs
  #TODO if pathway_ids are provided, get compound IDs that go with them
  # keggGetCompounds(pathway_id)
  # Download mols
  #TODO purrr is not already a dependency.  Try to convert this to lapply
  .get_mol_kegg <- function(compound_id) {
    mol <- KEGGREST::keggGet(compound_id, option = "mol")
    
    # Adds title to mol file because it is used later on by get_fx_groups()
    # currently get_fx_groups() queries KEGG for the title, but need to wean it of KEGG API
    names <- KEGGREST::keggGet(compound_id)[[1]]$NAME
    #just use the first name and remove separator
    title <- stringr::str_remove(names[1], ";")
    #add title line to mol file
    mol_clean <- gsub(">.*", "", mol) %>%
      paste0(title, "\n\n", .)
    mol_clean
  }
  mols <- lapply(compound_ids, .get_mol_kegg)

  file_paths <- fs::path(dir, compound_ids, ext = "mol")
  
  .write_mol <- function(mol_clean, file_path) {
    utils::write.table(
      mol_clean,
      file = file_path,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
  
  mapply(.write_mol, mol_clean = mols, file_path = file_paths)
  
  #construct tibble for output
  tibble::tibble(compound_id = compound_ids, mol_path = file_paths)
}


#' Get list of KEGG compound IDs for given KEGG pathway
#'
#' This is a temporary helper function until this function is improved and
#' pushed into KEGGREST package
#'
#' @param pathway string that is a KEGG identifier for a molecular pathway
keggGetCompounds <- function(pathway){
  
  resp <- 
    httr2::request("https://rest.kegg.jp/")  %>%  
    httr2::req_url_path("link/cpd/") %>%  
    httr2::req_url_path_append(pathway) %>%  
    httr2::req_perform()
  
  out <- resp %>% 
    httr2::resp_body_string() %>%  
    stringr::str_split_1("\n") %>%  
    stringr::str_extract("(?<=cpd:).*")
  out[!is.na(out)]

}