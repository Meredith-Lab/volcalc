#' Calculate volatility estimates for all compounds in a pathway
#'
#' Volatility value and category is estimated for all compounds in a specified
#' pathway using the SIMPOL formula
#'
#' @param pathway_id An optional character string specifying KEGG pathway ID, in
#'   format of 5 digits prepended with "map".
#' @param path An optional parameter to set relative path to location to
#'   download data.
#' @param redownload Download file again even if it has already been downloaded
#'   at path.
#' @param return_fx_groups When `TRUE`, includes functional group counts in
#'   final dataframe.
#' @param return_calc_steps When `TRUE`, includes intermediate volatility
#'   calculation steps in final dataframe.
#'
#' @return Dataframe with columns of basic compound info and volatility value
#'   and category. See documentation for column descriptions.
#'
#' @examples
#' \dontrun{
#' ex_pathway <- calc_pathway_vol(pathway_id = "map00361")
#' }
#' @export
calc_pathway_vol <-
  function(pathway_id,
           path = "data",
           redownload = FALSE,
           return_fx_groups = FALSE,
           return_calc_steps = FALSE) {
    
  compounds_from_pathway <- keggGetCompounds(pathway_id)
  mapply(save_compound_mol, compound_id = compounds_from_pathway,
         pathway_id = pathway_id, path = path, redownload = redownload)
  compound_files <- fs::dir_ls(fs::path(path, pathway_id))
  compound_names <- 
   fs::path_file(compound_files) %>% fs::path_ext_remove()
  
  #pre-allocate vectors for for-loop
  compounds_volatility <- c()
  compounds_fx_groups <- data.frame()
  for(compound in compound_names){
    compound_fx_groups <- get_fx_groups(compound_id = compound,
                                        pathway_id = pathway_id, path = path)
    compound_volatility <- calc_vol(compound_id = compound,
                                    pathway_id = pathway_id, path = path,
                                    save_file = FALSE, get_groups = FALSE,
                                    fx_groups_df = compound_fx_groups,
                                    return_fx_groups = return_fx_groups,
                                    return_calc_steps = return_calc_steps)
    compounds_volatility <- rbind(compounds_volatility, compound_volatility)
  }
  return(compounds_volatility)
}


