utils::globalVariables(".data")
#' SIMPOL.1 method for calculating estimated volatility
#' 
#' Implements the SIMPOL.1 group contribution method for predicting liquid vapor
#' pressure of organic compounds as described in Pankow & Asher (2008) and a
#' modified version described in Meredith et al. (2023).  Users will not usually
#' use this function directly, but rather through [calc_vol()] which uses this
#' as the default (currently only) method.
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
#' @param meredith logical; `FALSE`: use the original SIMPOL.1 method. `TRUE`:
#'   use the modified version in Meredith et al. (2023).
#'
#' @return the `fx_groups` tibble with the additional `log10_P` column 
#'
#' @references Pankow, J.F., Asher, W.E., 2008. SIMPOL.1: a simple group
#'   contribution method for predicting vapor pressures and enthalpies of
#'   vaporization of multifunctional organic compounds. Atmos. Chem. Phys.
#'   <https://doi.org/10.5194/acp-8-2773-2008>
#'   
#'   Meredith, L.K., Riemer, K., Geffre, P., Honeker, L., Krechmer, J., Graves,
#'   K., Tfaily, M., and Ledford, S.K. *In review*. Automating methods for
#'   estimating metabolite volatility.  Frontiers in Microbiology.
#'   
#' @seealso [calc_vol()]
#' @export
#'
#' @examples
#' mol_path <- mol_example("C16181.mol")
#' sdf <- ChemmineR::read.SDFset(mol_path)
#' fx_groups <- get_fx_groups(sdf)
#' simpol1(fx_groups)
simpol1 <- function(fx_groups, meredith = TRUE) {
  betas <- fx_groups %>%
    # assume NAs are 0s for the sake of this calculation
    dplyr::mutate(dplyr::across(dplyr::where(is.integer), function(x)
      ifelse(is.na(x), 0L, x))) %>%
    # TODO: There is also a table of coefs and equation to calculate b_k(T) at
    # different values of T. Why wasn't this implemented?
    dplyr::mutate(
      # multiplier for each functional group is volatility contribution to log10P
      # b_k(T)  * v_k,i
      b_00 = (1.79     * 1), #b_0(T) is an intercept/constant.  
      b_01 = (-0.438   * .data$carbons),
      b_02 = (-0.0338  * .data$carbons_asa),
      b_03 = (-0.675	 * .data$rings_aromatic),
      b_04 = (-0.0104  * .data$rings),
      b_05 = (-0.105   * .data$carbon_dbl_bonds),
      b_06 = (-0.506   * .data$CCCO_aliphatic_ring),
      b_07 = (-2.23	   * .data$hydroxyl_groups),
      b_08 = (-1.35	   * .data$aldehydes),
      b_09 = (-0.935   * .data$ketones),
      b_10 = (-3.58	   * .data$carbox_acids),
      b_11 = (-1.20	   * .data$ester),
      b_12 = (-0.718   * .data$ether),
      b_13 = (-0.683   * .data$ether_alicyclic),
      b_14 = (-1.03	   * .data$ether_aromatic),
      b_15 = (-2.23	   * .data$nitrate),
      b_16 = (-2.15	   * .data$nitro),
      b_17 = (-2.14	   * .data$hydroxyl_aromatic),
      b_18 = (-1.03	   * .data$amine_primary),
      b_19 = (-0.849	 * .data$amine_secondary),
      b_20 = (-0.608 	 * .data$amine_tertiary),
      b_21 = (-1.61    * .data$amine_aromatic),
      # TODO add amides to get_fx_groups() and uncomment these
      # b_22 = (-4.49    * .data$amide_primary),
      # b_23 = (-5.26    * .data$amide_secondary),
      # b_24 = (-2.63    * .data$amide_tertiary),
      b_25 = (-2.34    * .data$carbonylperoxynitrate),
      b_26 = (-0.368   * .data$peroxide),
      b_27 = (-2.48	   * .data$hydroperoxide),
      b_28 = (-2.48    * .data$carbonylperoxyacid),
      b_29 = (0.0432 	 * .data$nitrophenol),
      b_30 = (-2.67	   * .data$nitroester)
    ) 
  
  if (isTRUE(meredith)) {
    betas <- betas %>% 
      dplyr::mutate(
        # # Below are additions from Meredith et al.
        b_31 = (-2.23	  * .data$amides), #TODO remove once 1º, 2º, 3º amides are added
        b_32 = (-2.23	  * .data$phosphoric_acid),
        b_33 = (-2.23	  * .data$phosphoric_ester),
        b_34 = (-2.23	  * .data$sulfate),
        b_35 = (-2.23	  * .data$sulfonate),
        b_36 = (-2.23	  * .data$thiol),
        b_37 = (-1.20	  * .data$carbothioester)
      )
  }
  
  betas %>% 
    dplyr::mutate(log10_P = sum(dplyr::c_across(dplyr::starts_with("b_")))) %>% 
    dplyr::select(-dplyr::starts_with("b_"))
}