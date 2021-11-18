#' Calculate volatility estimate for compound
#'
#' Using functional group counts from get_fx_groups(), volatility is estimated
#' using the SIMPOL formula
#'
#' @param fx_groups_df dataframe of functional group counts for compounds
#'
#' @return input dataframe with new columns for volatility value and category
calc_vol <- function(fx_groups_df){
  aldehydes <- amine_aromatic <- amine_primary <- amine_secondary <- amine_tertiary <- carbon_dbl_bonds <- carbons <- carbox_acids <- case_when <- ester <- ether_alicyclic <- ether_aromatic <- hydroperoxide <- hydroxyl_groups <- ketones <- log_Sum <- log_alpha <- mass <- mutate <- nitrate <- nitro <- nitroester <- nitrophenol <- peroxide <- phenol <- rings <- rings_aromatic <- NULL
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
             (-0.94	 * ketones) %+%
             (-1.35	 * aldehydes) %+%
             (-2.23	 * hydroxyl_groups) %+%
             (-3.58	 * carbox_acids) %+%
             (-0.368 * peroxide) %+%
             (-2.48	 * hydroperoxide) %+%
             (-2.23	 * nitrate) %+%
             (-2.15	 * nitro) %+%
             (-0.1 	 * carbon_dbl_bonds) %+%
             (-0.01	 * rings) %+%
             (-0.68	 * rings_aromatic) %+%
             (-2.14	 * phenol) %+%
             (0.04 	 * nitrophenol) %+%
             (-2.67	 * nitroester) %+%
             (-1.2 	 * ester) %+%
             (-0.683 * ether_alicyclic) %+%
             (-1.03	 * ether_aromatic) %+%
             (-1.01	 * amine_primary) %+%
             (-0.85	 * amine_secondary) %+%
             (-0.6 	 * amine_tertiary) %+%
             (-1.61  * amine_aromatic),
           log_c = log_alpha + constant + log_Sum,
           volatility = dplyr::case_when(log_c <= 0 ~ "none",
                                  log_c > 0 & log_c <= 2 ~ "moderate",
                                  log_c > 2 ~ "high"))
  return(vol_df)
}
