#' Calculate volatility estimates for all compounds in a pathway
#'
#' Downloads all compound .mol files associated with a given KEGG pathway and
#' counts functional groups and calculates volatility value and category for all
#'
#' @param pathway_id character string that is 5 digits prepended with "map"
#' @param path optional parameter to set relative path to location to download data with default of creating "data" folder in home directory
#' @param redownload download files again even if they have already been downloaded at path
#' @param return_fx_groups whether to include columns of functional groups in final dataframe
#' @param return_calc_steps whether to include columns from intermediate volatility calculation steps in final dataframe
#'
#' @return dataframe with one row per compound and columns for basic compound info, functional group counts, and volatility estimate and category
#'
#' @export
calc_pathway_vol <- function(pathway_id, path = "data", redownload = FALSE,
                             return_fx_groups = FALSE, return_calc_steps = FALSE){
  compounds_volatility <- c()
  compounds_from_pathway <- keggGetCompounds(pathway_id)
  mapply(save_compound_mol, compound_id = compounds_from_pathway,
         pathway_id = pathway_id, path = path, redownload = redownload)
  compound_files <- list.files(paste0(path, "/", pathway_id))
  compound_names <- sapply(compound_files, function(x) substr(x, 1, nchar(x) - 4))
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

#' Get list of KEGG compound IDs for given KEGG pathway
#'
#' This is a temporary helper function until this function is improved and pushed into KEGGREST package
#'
#' @param pathway string that is a KEGG identifier for a molecular pathway
keggGetCompounds <- function(pathway){
  url <- sprintf("%s/link/cpd/%s", "http://rest.kegg.jp", pathway)
  .getUrl <- function (url, parser, ...)
  {
    url <- .cleanUrl(url)
    debug <- getOption("KEGGREST_DEBUG", FALSE)
    if (debug)
      .printf("url == %s", url)
    response <- httr::GET(url)
    httr::stop_for_status(response)
    content <- .strip(httr::content(response, "text"))
    if (nchar(content) == 0)
      return(character(0))
    do.call(parser, list(content, ...))
  }
  .printf <- function(...) message(noquote(sprintf(...)))
  .cleanUrl <- function (url)
  {
    url <- gsub(" ", "%20", url, fixed = TRUE)
    url <- gsub("#", "%23", url, fixed = TRUE)
    url <- gsub(":", "%3a", url, fixed = TRUE)
    sub("http(s)*%3a//", "http\\1://", url)
  }
  .strip <- function (str)
  {
    gsub("^\\s+|\\s+$", "", str)
  }
  .compoundParser <- function (txt)
  {
    cmptxt <- unlist(txt)
    lines <- strsplit(cmptxt, "\n")
    cmps <- gsub(".*cpd:", "", unlist(lines))
    cmps
  }
  .getUrl(url, .compoundParser)
}
