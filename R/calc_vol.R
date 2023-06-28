utils::globalVariables(".data")
#' Calculate volatility estimate for compound
#'
#' Volatility value and category is estimated for specified compound using the
#' SIMPOL formula
#'
#' @param input a path to a .mol file or a character representation such as SMILES or InChI 
#' @param from the form of `input` (currently only a path to a .mol file is implemented)
#' @param return_fx_groups When `TRUE`, includes functional group counts in
#'   final dataframe.
#' @param return_calc_steps When `TRUE`, includes intermediate volatility
#'   calculation steps in final dataframe.
#'
#' @return Dataframe with columns of basic compound info and volatility value
#'   and category. See documentation for column descriptions.
#'
#' @export
calc_vol <-
  function(input, 
           from = c("mol_path"),
           return_fx_groups = FALSE,
           return_calc_steps = FALSE) {
    
  from <- match.arg(from)
  
  #TODO: vectorize so `input` can be a vector
  if(from == "mol_path") {
    compound_sdf <- ChemmineR::read.SDFset(input)
  }
  fx_groups_df <- get_fx_groups(compound_sdf)
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
  
  vol_df %>% 
    dplyr::select(dplyr::all_of(c("formula", "name", "volatility", "category", cols_fx, cols_calc)))
  
  }

