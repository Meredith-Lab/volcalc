#' Calculate volatility estimate for compound
#'
#' Using functional group counts from get_fx_groups(), volatility is estimated
#' using the SIMPOL formula
#'
#' @param compound_id character string that is 5 digits prepended with a "C"
#' @param compound_formula character string of compound formula
#' @param pathway_id optional character string specifying KEGG pathway ID, in format of 5 digits prepended with "map"
#' @param path optional parameter to set relative path to location to download data with default of creating "data" folder in home directory
#' @param save_file save compound mol file using save_compound_mol function
#' @param redownload download file again even if it has already been downloaded at path
#' @param get_groups get dataframe of compound functional groups using get_fx_groups function
#' @param fx_groups_df dataframe of functional group counts for compounds, optional if reading functional groups dataframe in directly
#' @param return_fx_groups whether to include columns of functional groups in final dataframe
#' @param return_calc_steps whether to include columns from intermediate volatility calculation steps in final dataframe
#'
#' @return input dataframe with new columns for volatility value and category
#' @export
calc_vol <- function(compound_id = NULL, compound_formula = NULL, pathway_id = NULL,
                     path = "data",  redownload = FALSE, save_file = TRUE,
                     get_groups = TRUE, fx_groups_df = NULL,
                     return_fx_groups = FALSE, return_calc_steps = FALSE){
  if (is.null(compound_id) & is.null(compound_formula)) {
    stop("either compound_id or compound_formula needs to be specified")
  }
  if ((!isTRUE(save_file) & is.null(fx_groups_df)) | (!isTRUE(get_groups) & is.null(fx_groups_df)) |
      (!isTRUE(save_file) & !isTRUE(get_groups) & is.null(fx_groups_df))) {
    stop("either read in functional groups dataframe or set save_file and get_groups to true")
  }
  if (isTRUE(save_file)) {
    save_compound_mol(compound_id, compound_formula, pathway_id, path, redownload)
  }
  if (isTRUE(get_groups)) {
    fx_groups_df <- get_fx_groups(compound_id, pathway_id, path)
  }
  aldehydes <- amine_aromatic <- amine_primary <- amine_secondary <- amine_tertiary <- carbon_dbl_bonds <- carbons <- carbox_acids <- case_when <- ester <- ether_alicyclic <- ether_aromatic <- hydroperoxide <- hydroxyl_groups <- ketones <- log_Sum <- log_alpha <- mass <- mutate <- nitrate <- nitro <- nitroester <- nitrophenol <- peroxide <- phenol <- rings <- rings_aromatic <- amines <- amides <- phosphoric_acid <- phosphoric_ester <- sulfate <- sulfonate <- thiol <- carbothioester <- pathway <- name <- volatility <- category <- fluorines <- NULL
  # `constant` is vapor pressure baseline modified by functional group multipliers
  constant <- 1.79
  `%+%` <- function(x, y)  mapply(sum, x, y, MoreArgs = list(na.rm = TRUE))
  vol_df <- fx_groups_df %>%
    # mass is converted from grams to micrograms
    # 0.0000821 is universal gas constant
    # 293 is temperature in Kelvins
    dplyr::mutate(log_alpha = log((1000000*mass)/(0.0000821*293), base = 10),
           # multiplier for each functional group is volatility contribution
           log_Sum =
             (-0.44  * carbons) %+%
             (-0.935 * ketones) %+%
             (-1.35	 * aldehydes) %+%
             (-2.23	 * hydroxyl_groups) %+%
             (-3.58	 * carbox_acids) %+%
             (-0.368 * peroxide) %+%
             (-2.48	 * hydroperoxide) %+%
             (-2.23	 * nitrate) %+%
             (-2.15	 * nitro) %+%
             (-0.105 	 * carbon_dbl_bonds) %+%
             (-0.0104	 * rings) %+%
             (-0.675	 * rings_aromatic) %+%
             (-2.14	 * phenol) %+%
             (0.0432 	 * nitrophenol) %+%
             (-2.67	 * nitroester) %+%
             (-1.2 	 * ester) %+%
             (-0.683 * ether_alicyclic) %+%
             (-1.03	 * ether_aromatic) %+%
             (-1.03	 * amine_primary) %+%
             (-0.849	 * amine_secondary) %+%
             (-0.608 	 * amine_tertiary) %+%
             (-1.61  * amine_aromatic) %+%
             (-2.23	 * amines) %+%
             (-2.23	 * amides) %+%
             (-2.23	 * phosphoric_acid) %+%
             (-2.23	 * phosphoric_ester) %+%
             (-2.23	 * sulfate) %+%
             (-2.23	 * sulfonate) %+%
             (-2.23	 * thiol) %+%
             (-1.2 	 * carbothioester),
           volatility = log_alpha + constant + log_Sum,
           category = dplyr::case_when(volatility <= 0 ~ "low",
                                  volatility > 0 & volatility <= 2 ~ "moderate",
                                  volatility > 2 ~ "high"))
  if (isTRUE(return_fx_groups) & !isTRUE(return_calc_steps)){
    subset_vol_df <- dplyr::select(vol_df, pathway:name, volatility:category, carbons:fluorines)
  } else if (!isTRUE(return_fx_groups) & isTRUE(return_calc_steps)){
    subset_vol_df <- dplyr::select(vol_df, pathway:name, volatility:category, mass, log_alpha:log_Sum)
  } else if (isTRUE(return_fx_groups) & isTRUE(return_calc_steps)) {
    subset_vol_df <- vol_df
  } else {
    subset_vol_df <- dplyr::select(vol_df, pathway:name, volatility:category)
  }
  return(subset_vol_df)
}
