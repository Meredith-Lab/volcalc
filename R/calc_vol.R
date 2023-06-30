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
#' @return a tibble with columns of basic compound info and volatility value
#'   and category. See documentation for column descriptions.
#'
#' @export
calc_vol <-
  function(input, 
           from = c("mol_path"),
           return_fx_groups = FALSE,
           return_calc_steps = FALSE) {
  #TODO: consider bringing back the fx_groups_df argument to pass a manually constructed dataframe since not all functional groups in Pankow & Asher are detected by get_fx_groups currently.
    
  from <- match.arg(from)
  
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

