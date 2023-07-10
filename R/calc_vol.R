#' Calculate volatility estimate for a compound
#'
#' Relative volatility value and category is estimated for specified compound
#' using group contribution methods.
#'
#' @param input a path to a .mol file or a character representation such as
#'   SMILES or InChI
#' @param from the form of `input` (currently only a path to a .mol file is
#'   implemented)
#' @param method the method for calculating estimated volatility. Currently only
#'   the SIMPOL.1 method is implemented---see [simpol1()] for more details.
#' @param return_fx_groups When `TRUE`, includes functional group counts in
#'   final dataframe.
#' @param return_calc_steps When `TRUE`, includes intermediate volatility
#'   calculation steps in final dataframe.
#' 
#'
#' @return a tibble with columns of basic compound info and volatility value
#'   and category created by [get_fx_groups()] and [simpol1()].
#'   
#' @seealso [get_fx_groups()], [simpol1()]
#'
#' @export
calc_vol <-
  function(input, 
           from = c("mol_path"),
           method = c("simpol1"),
           return_fx_groups = FALSE,
           return_calc_steps = FALSE) {
    
  from <- match.arg(from)
  #for future extensions in case other methods are added
  method <- match.arg(method)
  
  if(from == "mol_path") {
    compound_sdf_list <- lapply(input, ChemmineR::read.SDFset)
  }
  
  fx_groups_df <- lapply(compound_sdf_list, get_fx_groups) %>% 
    dplyr::bind_rows()
  
  # calculate volatility
  vol_df <- simpol1(fx_groups_df) 
  
  # wrangle output
  cols_fx <- NULL
  cols_calc <- NULL
  if (isTRUE(return_fx_groups)) {
    cols_fx <- colnames(fx_groups_df)[!colnames(fx_groups_df) %in% c("formula", "name", "mass")]
  }
  if (isTRUE(return_calc_steps)) {
    cols_calc <- c("mass", "log_alpha", "log_Sum")
  }
  
  #return:
  vol_df %>% 
    dplyr::select(dplyr::all_of(c("formula", "name", "volatility", "category", cols_fx, cols_calc)))
  }

