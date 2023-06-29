utils::globalVariables(".data")
#' SIMPOL.1 method for calculating estimated volatility
#' 
#' Implements SIMPOL.1
#'
#' @param fx_groups a data.frame with counts of functional groups produced by
#'   `get_fx_groups()`
#'
#' @return a `data.frame` object
#' @export
#'
#' @examples
simpol1 <- function(fx_groups) {
  
  # `constant` is vapor pressure baseline modified by functional group multipliers
  constant <- 1.79 # b_0(T)
  vol_df <- 
    fx_groups %>%
    # mass is converted from grams to micrograms
    # 0.0000821 is universal gas constant
    # 293 is temperature in Kelvins (19.85ºC)
    # TODO: not sure why hard-coded at 293K.  The coefficients below are from a
    # table that used 293.15K (table 6 in Pankow & Asher 2008).  There is also a
    # table of coefs and equation to calculate b_k(T) at different values of T.
    # Why wasn't this implemented?
    dplyr::mutate(
      log_alpha = log((1000000 * .data$mass) / (0.0000821 * 293), base = 10),
      # multiplier for each functional group is volatility contribution
      log_Sum = #TODO why is this called log_Sum?  There's no log operation, right?
        #TODO deal with NAs before this step and get rid of %+% operator
        # b_k(T)  * v_k,i
        (-0.438   * .data$carbons) %+%
        (-0.935   * .data$ketones) %+%
        (-1.35	  * .data$aldehydes) %+%
        (-2.23	  * .data$hydroxyl_groups) %+%
        (-3.58	  * .data$carbox_acids) %+%
        (-0.368   * .data$peroxide) %+%
        (-2.48	  * .data$hydroperoxide) %+%
        (-2.23	  * .data$nitrate) %+%
        (-2.15	  * .data$nitro) %+%
        (-0.105   * .data$carbon_dbl_bonds) %+%
        (-0.0104  * .data$rings) %+%
        (-0.675	  * .data$rings_aromatic) %+%
        (-2.14	  * .data$phenol) %+%
        (0.0432 	* .data$nitrophenol) %+%
        (-2.67	  * .data$nitroester) %+%
        (-1.20	  * .data$ester) %+%
        (-0.718   * .data$ether) %+%
        (-0.683   * .data$ether_alicyclic) %+%
        (-1.03	  * .data$ether_aromatic) %+%
        (-1.03	  * .data$amine_primary) %+%
        (-0.849	  * .data$amine_secondary) %+%
        (-0.608 	* .data$amine_tertiary) %+%
        (-1.61    * .data$amine_aromatic) %+%
        (-2.23	  * .data$amines) %+%
        (-2.23	  * .data$amides) %+%
        (-2.23	  * .data$phosphoric_acid) %+%
        (-2.23	  * .data$phosphoric_ester) %+%
        (-2.23	  * .data$sulfate) %+%
        (-2.23	  * .data$sulfonate) %+%
        (-2.23	  * .data$thiol) %+%
        (-1.20	  * .data$carbothioester),
      
      #TODO shoud the following be part of simpol1() or part of calc_vol() ?
      # I think constant + log_Sum = log10Pº_{L,i}(T), not sure what adding log_alpha makes it
      volatility = constant + .data$log_alpha + .data$log_Sum, #units are atm ?? (or log10 atm maybe)
      category = dplyr::case_when(
        .data$volatility <  -2                        ~ "non",
        .data$volatility >= -2 & .data$volatility < 0 ~ "low",
        .data$volatility >= 0  & .data$volatility < 2 ~ "intermediate",
        .data$volatility >= 2                         ~ "high"
      )
    ) 
}