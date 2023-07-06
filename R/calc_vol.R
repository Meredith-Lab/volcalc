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
#' @param fx_groups optional data frame in the format produced by
#'   [get_fx_groups()].  If provided, it will override calculated functional
#'   group counts---use with caution and see Details for more.
#' @param return_fx_groups When `TRUE`, includes functional group counts in
#'   final dataframe.
#' @param return_calc_steps When `TRUE`, includes intermediate volatility
#'   calculation steps in final dataframe.
#'   
#' @details By default, [calc_vol()] uses [get_fx_groups()] to get counts of
#' functional groups and other properties that contribute to volatility.
#' However, not all functional groups are currently implemented, and it's
#' possible that functional group counts may not be captured perfectly for all
#' molecules.  Therefore, you may provide a data frame of manual counts to
#' `fx_groups`.  It expects one row per compound, and exactly the column names
#' produced by [get_fx_groups()].
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
           fx_groups = NULL,
           return_fx_groups = FALSE,
           return_calc_steps = FALSE) {
    
  from <- match.arg(from)
  #for future extensions in case other methods are added
  method <- match.arg(method)
  
  if(from == "mol_path") {
    compound_sdf_list <- lapply(input, ChemmineR::read.SDFset)
  }
  if (is.null(fx_groups)) {
    fx_groups_df <- lapply(compound_sdf_list, get_fx_groups) %>% 
      dplyr::bind_rows()
  } else {
    #TODO validate fx_groups
    #TODO write tests for this argument
    fx_groups_df <- fx_groups
  }
  
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

