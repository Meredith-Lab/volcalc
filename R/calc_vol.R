#' Calculate volatility estimate for compound
#'
#' Volatility value and category is estimated for specified compound using the
#' SIMPOL formula
#'
#' @param compound_id A character string that is 5 digits prepended with a "C".
#' @param compound_formula A character string detailing a compound formula.
#' @param pathway_id An optional character string specifying KEGG pathway ID, in
#'   format of 5 digits prepended with "map".
#' @param path An optional parameter to set relative path to location to
#'   download data.
#' @param save_file Whether to save downloaded compound mol files.
#' @param redownload Download file again even if it has already been downloaded
#'   at path.
#' @param get_groups When `FALSE`, will expect a dataframe to be read in with
#'   `fx_groups_df` argument.
#' @param fx_groups_df A dataframe of functional group counts for compounds
#'   generated from [`get_fx_groups()`].
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
#' ex_compound <- calc_vol(compound_id = "C16181")
#' }
#' @export
calc_vol <-
  function(compound_id = NULL,
           compound_formula = NULL,
           pathway_id = NULL,
           path = "data",
           redownload = FALSE,
           save_file = TRUE,
           get_groups = TRUE,
           fx_groups_df = NULL,
           return_fx_groups = FALSE,
           return_calc_steps = FALSE) {
    
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
    
  #assign variables to quiet devtools::check()
  aldehydes <- amine_aromatic <- amine_primary <- amine_secondary <- amine_tertiary <- carbon_dbl_bonds <- carbons <- carbox_acids <- case_when <- ester <- ether <- ether_alicyclic <- ether_aromatic <- hydroperoxide <- hydroxyl_groups <- ketones <- log_Sum <- log_alpha <- mass <- mutate <- nitrate <- nitro <- nitroester <- nitrophenol <- peroxide <- hydroxyl_aromatic <- rings <- rings_aromatic <- amides <- phosphoric_acid <- phosphoric_ester <- sulfate <- sulfonate <- thiol <- carbothioester <- pathway <- name <- volatility <- category <- fluorines <- NULL
  
  # `constant` is vapor pressure baseline modified by functional group multipliers
  constant <- 1.79
  vol_df <- fx_groups_df %>%
    # mass is converted from grams to micrograms
    # 0.0000821 is universal gas constant
    # 293 is temperature in Kelvins
    dplyr::mutate(log_alpha = log((1000000*mass)/(0.0000821*293.15), base = 10),
           # multiplier for each functional group is volatility contribution
           log_Sum =
             (-0.438  * carbons) %+%
             (-0.935  * ketones) %+%
             (-1.35	  * aldehydes) %+%
             (-2.23	  * hydroxyl_groups) %+%
             (-3.58	  * carbox_acids) %+%
             (-0.368  * peroxide) %+%
             (-2.48	  * hydroperoxide) %+%
             (-2.23	  * nitrate) %+%
             (-2.15	  * nitro) %+%
             (-0.105  * carbon_dbl_bonds) %+%
             (-0.0104 * rings) %+%
             (-0.675	* rings_aromatic) %+%
             (-2.14	  * hydroxyl_aromatic) %+%
             (0.0432 	* nitrophenol) %+%
             (-2.67	  * nitroester) %+%
             (-1.20	  * ester) %+%
             (-0.718  * ether) %+%
             (-0.683  * ether_alicyclic) %+%
             (-1.03	  * ether_aromatic) %+%
             (-1.03	  * amine_primary) %+%
             (-0.849	* amine_secondary) %+%
             (-0.608 	* amine_tertiary) %+%
             (-1.61   * amine_aromatic) %+%
             (-2.23	  * amides) %+%
             (-2.23	  * phosphoric_acid) %+%
             (-2.23	  * phosphoric_ester) %+%
             (-2.23	  * sulfate) %+%
             (-2.23	  * sulfonate) %+%
             (-2.23	  * thiol) %+%
             (-1.20	  * carbothioester),
           volatility = log_alpha + constant + log_Sum,
           category = dplyr::case_when(
             volatility <  -2                  ~ "non", #less than -2
             volatility >= -2 & volatility < 0 ~ "low", #great than or equal to -2, less than 0
             volatility >= 0 & volatility < 2  ~ "intermediate", #greater than or equal to 0, less than 2
             volatility >= 2                   ~ "high" #greater than or equal to 2
           )) 
  if (isTRUE(return_fx_groups) & !isTRUE(return_calc_steps)){
    subset_vol_df <-
      dplyr::select(vol_df, pathway:name, volatility:category, carbons:fluorines)
  } else if (!isTRUE(return_fx_groups) & isTRUE(return_calc_steps)){
    subset_vol_df <- 
      dplyr::select(vol_df, pathway:name, volatility:category, mass, log_alpha:log_Sum)
  } else if (isTRUE(return_fx_groups) & isTRUE(return_calc_steps)) {
    subset_vol_df <- vol_df
  } else {
    subset_vol_df <- dplyr::select(vol_df, pathway:name, volatility:category)
  }
  return(subset_vol_df)
}
