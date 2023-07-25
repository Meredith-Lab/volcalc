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
#' @return a tibble with relative volatility index (`volatility`) and volatility
#'   category (`category`).
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
  
  #TODO: needs testing before implementing
  # if(from == "smiles") { 
  #   compound_sdf_list <- lapply(input, ChemmineR::smiles2sdf)
  # }
  
  fx_groups_df <-
    lapply(compound_sdf_list, get_fx_groups) %>% 
    dplyr::bind_rows(.id = {{from}}) #adds column for input named "mol_path" or "smiles"
  
  # calculate relative volatility & categories from logP
  vol_df <- simpol1(fx_groups_df) %>% 
    dplyr::mutate(
      # mass is converted from grams to micrograms
      # 0.0000821 is universal gas constant
      # 293.15 is temperature in Kelvins (20ºC)
      log_alpha = log10((1000000 * .data$mass) / (0.0000821 * 293.15)),
      volatility = .data$log_alpha + .data$log10_P, 
      #TODO add @details to documentation explaining categories.  Make sure they match manuscript
      category = dplyr::case_when(
        .data$volatility <  -2                        ~ "non",
        .data$volatility >= -2 & .data$volatility < 0 ~ "low",
        .data$volatility >= 0  & .data$volatility < 2 ~ "intermediate",
        .data$volatility >= 2                         ~ "high"
      )
    )
  
  # wrangle output
  cols_fx <- NULL
  cols_calc <- NULL
  if (isTRUE(return_fx_groups)) {
    cols_fx <- colnames(fx_groups_df)[!colnames(fx_groups_df) %in% c("formula", "name", "mass")]
  }
  if (isTRUE(return_calc_steps)) {
    #TODO document log_alpha (need to figure out why it's called that first)
    cols_calc <- c("mass", "log_alpha", "log10_P")
  }
  
  #return:
  vol_df %>% 
    dplyr::select(dplyr::all_of(c({{from}}, "formula", "name", "volatility", "category", cols_fx, cols_calc)))

  }

