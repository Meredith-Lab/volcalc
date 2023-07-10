utils::globalVariables(".data")
#' SIMPOL.1 method for calculating estimated volatility
#' 
#' Implements the SIMPOL.1 group contribution method for predicting liquid vapor
#' pressure of organic compounds as described in Pankow & Asher (2008).  Users
#' will not usually use this function directly, but rather through [calc_vol()]
#' which uses this as the default (currently only) method.
#' 
#'
#' @param fx_groups a data.frame or tibble with counts of functional groups
#'   produced by [get_fx_groups()] (or manually, with the same column names)
#'
#' @return The `fx_groups` tibble with the additional columns:
#'   * `log_alpha` — not sure! units?
#'   * `log_Sum` — \eqn{\sum_k\nu_{k,i}b_k(T)}, or the sum of the counts of functional groups (\eqn{\nu_{k,i}}) times the coefficients for each functional group (\eqn{b_K(T)}). units?
#'   * `volatility` — relative volatility? units?
#'   * `category` — category of volatility?
#'
#' @references Pankow, J.F., Asher, W.E., 2008. SIMPOL.1: a simple group
#'   contribution method for predicting vapor pressures and enthalpies of
#'   vaporization of multifunctional organic compounds. Atmos. Chem. Phys.
#'   https://doi.org/10.5194/acp-8-2773-2008
#' @seealso [calc_vol()]
#' @export
#'
#' @examples
#' mol_path <- mol_example("C16181.mol")
#' sdf <- ChemmineR::read.SDFset(mol_path)
#' fx_groups <- get_fx_groups(sdf)
#' simpol1(fx_groups)
simpol1 <- function(fx_groups) {
  fx_groups %>%
    # assume NAs are 0s for the sake of this calculation
    dplyr::mutate(dplyr::across(dplyr::where(is.integer), function(x)
      ifelse(is.na(x), 0L, x))) %>%
    
    
    # TODO: There is also a table of coefs and equation to calculate b_k(T) at
    # different values of T. Why wasn't this implemented?
    dplyr::mutate(
      # multiplier for each functional group is volatility contribution
      log_Sum = #TODO why is this called log_Sum?  There's no log operation, right?
        #TODO I think log_Sum should include `constant`, as it is equivalent to b_0(T) * 1
        # b_k(T)  * v_k,i
        (1.79     * 1) + #b_0(T) is an intercept/constant
        (-0.438   * .data$carbons) +
        (-0.935   * .data$ketones) +
        (-1.35	  * .data$aldehydes) +
        (-2.23	  * .data$hydroxyl_groups) +
        (-3.58	  * .data$carbox_acids) +
        (-0.368   * .data$peroxide) +
        (-2.48	  * .data$hydroperoxide) +
        (-2.23	  * .data$nitrate) +
        (-2.15	  * .data$nitro) +
        (-0.105   * .data$carbon_dbl_bonds) +
        (-0.0104  * .data$rings) +
        (-0.675	  * .data$rings_aromatic) +
        (-2.14	  * .data$phenol) +
        (0.0432 	* .data$nitrophenol) +
        (-2.67	  * .data$nitroester) +
        (-1.20	  * .data$ester) +
        (-0.718   * .data$ether) +
        (-0.683   * .data$ether_alicyclic) +
        (-1.03	  * .data$ether_aromatic) +
        (-1.03	  * .data$amine_primary) +
        (-0.849	  * .data$amine_secondary) +
        (-0.608 	* .data$amine_tertiary) +
        (-1.61    * .data$amine_aromatic) +
        (-2.23	  * .data$amines) + #TODO doesn't this double count amines?
        (-2.23	  * .data$amides) +
        (-2.23	  * .data$phosphoric_acid) +
        (-2.23	  * .data$phosphoric_ester) +
        (-2.23	  * .data$sulfate) +
        (-2.23	  * .data$sulfonate) +
        (-2.23	  * .data$thiol) +
        (-1.20	  * .data$carbothioester)
    ) 
}