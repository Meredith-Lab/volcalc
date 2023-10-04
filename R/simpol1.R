utils::globalVariables(".data")
#' SIMPOL.1 method for calculating estimated volatility
#' 
#' Implements the SIMPOL.1 group contribution method for predicting liquid vapor
#' pressure of organic compounds as described in Pankow & Asher (2008).  Users
#' will not usually use this function directly, but rather through [calc_vol()]
#' which uses this as the default (currently only) method.
#' 
#' @details The output includes a column for `log10_P` where
#' \eqn{\textrm{log}_{10} P_{\textrm{L},i}^\circ(T) = \sum_k\nu_{k,i}b_k(T)}, or
#' the sum of the counts of functional groups (\eqn{\nu_{k,i}}) times the
#' coefficients for each functional group (\eqn{b_K(T)}). Units are in log10
#' atmospheres.
#' 
#' @note The method described in Pankow & Asher (2008) allows for
#' calculations of logP at different temperatures.  This implementation
#' currently only calculates values at 20ºC.
#' 
#' @param fx_groups a data.frame or tibble with counts of functional groups
#'   produced by [get_fx_groups()] (or manually, with the same column names)
#'
#' @return the `fx_groups` tibble with the additional `log10_P` column 
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
  # TODO: this is actually a *modified* version of SIMPOL.1 with additional
  # functional group coefs added.  Should there be a separate function for the
  # original SMIPOL.1 method that doesn't count these additional groups
  fx_groups %>%
    # assume NAs are 0s for the sake of this calculation
    dplyr::mutate(dplyr::across(dplyr::where(is.integer), function(x)
      ifelse(is.na(x), 0L, x))) %>%
    # TODO: There is also a table of coefs and equation to calculate b_k(T) at
    # different values of T. Why wasn't this implemented?
    dplyr::mutate(
      # multiplier for each functional group is volatility contribution to log10P
      log10_P = 
        # b_k(T)  * v_k,i
        (1.79     * 1) + #b_0(T) is an intercept/constant.  
        (-0.438   * .data$carbons) +
        (-0.0338  * .data$carbons_asa) +
        (-0.675	  * .data$rings_aromatic) +
        (-0.0104  * .data$rings) +
        (-0.105   * .data$carbon_dbl_bonds) +
        (-0.506   * .data$CCCO_aliphatic_ring) +
        (-2.23	  * .data$hydroxyl_aliphatic) +
        (-1.35	  * .data$aldehydes) +
        (-0.935   * .data$ketones) +
        (-3.58	  * .data$carbox_acids) +
        (-1.20	  * .data$ester) +
        (-0.718   * .data$ether) +
        (-0.683   * .data$ether_alicyclic) +
        (-1.03	  * .data$ether_aromatic) +
        (-2.23	  * .data$nitrate) +
        (-2.15	  * .data$nitro) +
        (-2.14	  * .data$hydroxyl_aromatic) +
        (-1.03	  * .data$amine_primary) +
        (-0.849	  * .data$amine_secondary) +
        (-0.608 	* .data$amine_tertiary) +
        (-1.61    * .data$amine_aromatic) +
        #TODO: missing amides
        (-2.34    * .data$carbonylperoxynitrate) +
        (-0.368   * .data$peroxide) +
        (-2.48	  * .data$hydroperoxide) +
        (-2.48    * .data$carbonylperoxyacid) +
        (0.0432 	* .data$nitrophenol) +
        (-2.67	  * .data$nitroester) +
        # Below are additions from Meredith et al.
        (-2.23	  * .data$amides) + #TODO remove once 1º, 2º, 3º amides are added
        (-2.23	  * .data$phosphoric_acid) +
        (-2.23	  * .data$phosphoric_ester) +
        (-2.23	  * .data$sulfate) +
        (-2.23	  * .data$sulfonate) +
        (-2.23	  * .data$thiol) +
        (-1.20	  * .data$carbothioester)
    ) 
}