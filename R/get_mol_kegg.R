utils::globalVariables(".data")
#' Download compound .mol files from KEGG
#'
#' Downloads mol files corresponding to individual compounds or compounds in a
#' pathway from KEGG.
#'
#' @param compound_ids character vector of KEGG compound IDs---5 digits
#'   prepended with a "C".
#' @param pathway_ids character vector of KEGG pathway or pathway module IDs---5
#'   digits prepended with "map" or "M", respectively.
#' @param dir path to a folder to save .mol files in. Folder will be created if
#'   it does not already exist
#' @param force logical; by default (`FALSE`), .mol files will not be downloaded
#'   if they are found in `dir`.  Set this to `TRUE` to download and overwrite
#'   existing files.
#'
#' @returns a tibble with the columns `compound_ids`, `pathway_ids` (if used),
#'   and `mol_paths` (paths to downloaded .mol files)
#' @export
#'
#' @examples
#' \dontrun{
#' get_mol_kegg(compound_ids = c("C16181", "C06074"), dir = tempdir())
#' get_mol_kegg(pathway_ids = "map00253", dir = tempdir())
#' }
get_mol_kegg <- function(compound_ids, pathway_ids, dir, force = FALSE){
  
  if(missing(dir)) stop("`dir` is required")
  if ((missing(compound_ids) & missing(pathway_ids)) |
      !missing(compound_ids) & !missing(pathway_ids)) {
    stop("One of `compound_id` or `pathway_id` are required")
  }
  
  #if compounds are provided
  if (!missing(compound_ids)) {
    if (!all(stringr::str_detect(compound_ids, "^[C][:digit:]{5}$"))) {
      stop("Some compound_ids are not in the correct KEGG format")
    }
    fs::dir_create(dir)
    out_tbl <-
      tibble::tibble(compound_id = compound_ids) %>% 
      dplyr::mutate(mol_path = fs::path(dir, .data$compound_id, ext = "mol"))
  }
  # if pathways are provided
  if (!missing(pathway_ids)) {
    if (!all(stringr::str_detect(pathway_ids, "^(map|M)\\d{5}$"))) {
      stop("Some pathway_ids are not in the correct KEGG format")
    }
    fs::dir_create(dir, pathway_ids)
    compound_ids_list <- lapply(pathway_ids, keggGetCompounds)
    names(compound_ids_list) <- pathway_ids
    out_tbl <- 
      tibble::enframe(compound_ids_list, name = "pathway_id", value = "compound_id") %>% 
      tidyr::unnest(tidyselect::everything()) %>% 
      dplyr::mutate(mol_path = fs::path(dir, .data$pathway_id, .data$compound_id, ext = "mol"))
    
  }
  
  if(isFALSE(force)) {
    to_dl <- out_tbl$compound_id[!fs::file_exists(out_tbl$mol_path)]
    out_paths <- out_tbl$mol_path[!fs::file_exists(out_tbl$mol_path)]
  } else {
    to_dl <- out_tbl$compound_id
    out_paths <- out_tbl$mol_path
  }
  
  if (length(to_dl) == 0) {
    #if nothing to download, return early
    return(out_tbl)
  } else {
    
    # Download mols
    .get_mol_kegg <- function(compound_id) {
      #TODO I think the KEGG API can handle up to 10 requests at once separated by
      #"+".  Unfortunately all the mol files come out in a single textfile and
      #would need parsing to separate.  Could speed things up by reducing API
      #calls, but would require additional code.
      mol <- KEGGREST::keggGet(compound_id, option = "mol")
      
      # Adds title to mol file because it is used later on by get_fx_groups()
      names <- KEGGREST::keggGet(compound_id)[[1]]$NAME
      # Only use the first name and remove separator
      title <- stringr::str_remove(names[1], ";")
      # add title line to mol file
      mol_clean <- paste0(title, "\n\n\n", gsub(">.*", "", mol))
      mol_clean
    }
    mols <- lapply(to_dl, .get_mol_kegg)
    
    # write mol files
    .write_mol <- function(mol_clean, file_path) {
      utils::write.table(
        mol_clean,
        file = file_path,
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
      )
    }
    
    mapply(.write_mol, mol_clean = mols, file_path = out_paths)
    
    return(out_tbl)
  }
}


#' Get list of KEGG compound IDs for given KEGG pathway
#'
#' This is a temporary helper function until this function is improved and
#' pushed into KEGGREST package
#'
#' @param pathway string that is a KEGG identifier for a molecular pathway
#' @noRd
keggGetCompounds <- function(pathway){
  
  resp <- 
    httr2::request("https://rest.kegg.jp/")  %>%  
    httr2::req_url_path("link/cpd/") %>%  
    httr2::req_url_path_append(pathway) %>% 
    httr2::req_retry(max_tries = 3) %>% 
    httr2::req_perform()
  
  out <- resp %>% 
    httr2::resp_body_string() %>%  
    stringr::str_split_1("\n") %>%  
    stringr::str_extract("(?<=cpd:).*")
  out[!is.na(out)]
  
}



