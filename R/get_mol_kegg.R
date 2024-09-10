utils::globalVariables(".data")
#' Download compound .mol files from KEGG
#'
#' Downloads mol files corresponding to individual compounds or compounds in a
#' pathway from KEGG.
#'
#' @param compound_ids Character vector of KEGG compound IDs—5 digits
#'   prepended with a "C".
#' @param pathway_ids Character vector of KEGG pathway or pathway module IDs—5
#'   digits prepended with "map" or "M", respectively.
#' @param dir Path to a folder to save .mol files in. Folder will be created if
#'   it does not already exist.
#' @param force Logical; by default (`FALSE`), .mol files will not be downloaded
#'   if they are found in `dir`. Set this to `TRUE` to download and overwrite
#'   existing files.
#'   
#' @note For additional functionality for interacting with KEGG, try the
#'   `KEGGREST` package, which this function was inspired by.
#'   
#' @returns A tibble with the columns `compound_ids`, `pathway_ids` (if used),
#'   and `mol_paths` (paths to downloaded .mol files).
#'   
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
    stop("One of `compound_ids` or `pathway_ids` are required")
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
    compound_ids_list <- lapply(pathway_ids, get_compounds_kegg)
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
    mols <- dl_mol_kegg(to_dl)
    
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
#' @param pathway string that is a KEGG identifier for a molecular pathway
#' @noRd
get_compounds_kegg <- function(pathway){
  
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

#' Get and wrangle mol files for a single API request of up to 10 IDs
#' @noRd
.dl_mol_kegg <- function(ids) {
  if (length(ids) > 10) {
    stop("Provide 10 or fewer IDs at a time")
  }
  req_names <- 
    httr2::request("https://rest.kegg.jp/get") %>% 
    httr2::req_url_path_append(paste(ids, collapse = "+")) %>% 
    httr2::req_retry(max_tries = 3)

  resp_names <- httr2::req_perform(req_names) %>% 
    httr2::resp_body_string()
  
  # There's a lot of stuff in the response, but I only care about the compound name
  names <- resp_names %>%
    stringr::str_extract_all("(?<=NAME).+(?=\\n)") %>%
    unlist() %>% 
    stringr::str_trim() %>% 
    stringr::str_remove(";")
  
  # get mol file
  req_mols <- req_names %>% 
    httr2::req_url_path_append("mol")
  
  resp_mols <- httr2::req_perform(req_mols) %>% 
    httr2::resp_body_string()
  
  # wrangle into valid mol files
  mols <- resp_mols %>% 
    stringr::str_split("(?<=\\${4})", n = length(ids)) %>%
    unlist() %>% 
    stringr::str_trim(side = "left")
  mols <- 
    gsub(">.*", "", mols) #for some reason this pattern doesn't work with str_remove()
  
  #add compound name in correct place
  paste0(names, "\n\n\n", mols)
}

#' Split vector into list elements of max length
#' @noRd
split_to_list <- function(x, max_len = 10) {
  
  if(length(x) > max_len) {
    n_groups <- ceiling(length(x) / max_len)
    split(x, f = cut(seq_along(x), breaks = n_groups)) %>%
      purrr::set_names(NULL)
  } else {
    list(x)
  }
  
}

#' Get mol files for compound_ids by splitting into groups of 10 and calling .dl_mol_kegg
#' @noRd
dl_mol_kegg <- function(compound_ids) {
  #balances compound_ids into groups of less than 10 to meet API guidelines
  compound_id_list <- split_to_list(compound_ids, max_len = 10)
  
  #maps over list, but returns it to a single character vector to simplify wrangling code
  purrr::map(compound_id_list, .dl_mol_kegg) %>% 
    purrr::list_c() %>% 
    glue::glue_collapse()
}
