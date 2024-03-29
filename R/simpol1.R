#' SIMPOL.1 method for calculating estimated volatility
#' 
#' Implements the SIMPOL.1 group contribution method for predicting liquid vapor
#' pressure of organic compounds as described in Pankow & Asher (2008) and a
#' modified version described in Meredith et al. (2023).  Users will not usually
#' use this function directly, but rather through [calc_vol()].
#' 
#' @details The output includes a column for `log10_P` where
#' \eqn{\textrm{log}_{10} P_{\textrm{L},i}^\circ(T) = \sum_k\nu_{k,i}b_k(T)}, or
#' the sum of the counts of functional groups (\eqn{\nu_{k,i}}) times the
#' coefficients for each functional group (\eqn{b_K(T)}). Units are in log10
#' atmospheres.
#' 
#' The modified method in Meredith et al. (2023) adds the following additional
#' functional groups and coefficients:
#' 
#' - Phosphoric acid (-2.23)
#' - Phosphoric ester (-2.23)
#' - Sulfate (-2.23)
#' - Sulfonate (-2.23)
#' - Thiol (-2.23)
#' - Carbothioester (-1.20)
#' 
#' @note The method described in Pankow & Asher (2008) allows for
#' calculations of logP at different temperatures.  This implementation
#' currently only calculates values at 20ºC.
#' 
#' @param fx_groups A data.frame or tibble with counts of functional groups
#'   produced by [get_fx_groups()] (or manually, with the same column names).
#' @param meredith Logical; `FALSE`: use the original SIMPOL.1 method. `TRUE`:
#'   use the modified version in Meredith et al. (2023).
#'
#' @returns The `fx_groups` tibble with the additional `log10_P` column.
#'
#' @references 
#'   Meredith L, Ledford S, Riemer K, Geffre P, Graves K, Honeker L, LeBauer D,
#'   Tfaily M, Krechmer J. 2023. Automating methods for estimating metabolite
#'   volatility. Frontiers in Microbiology. \doi{10.3389/fmicb.2023.1267234}
#'   
#'   Pankow, J.F., Asher, W.E. 2008. SIMPOL.1: a simple group
#'   contribution method for predicting vapor pressures and enthalpies of
#'   vaporization of multifunctional organic compounds. Atmos. Chem. Phys.
#'   \doi{10.5194/acp-8-2773-2008}
#'   
#' @seealso [calc_vol()]
#' @export
#'
#' @examples
#' mol_path <- mol_example()[3]
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
      b_04 = (-0.0104  * .data$rings_aliphatic),
      b_05 = (-0.105   * .data$carbon_dbl_bonds_aliphatic),
      b_06 = (-0.506   * .data$CCCO_aliphatic_ring),
      b_07 = (-2.23	   * .data$hydroxyl_aliphatic),
      b_08 = (-1.35	   * .data$aldehydes),
      b_09 = (-0.935   * .data$ketones),
      b_10 = (-3.58	   * .data$carbox_acids),
      b_11 = (-1.20	   * .data$ester),
      b_12 = (-0.718   * .data$ether_alkyl),
      b_13 = (-0.683   * .data$ether_alicyclic),
      b_14 = (-1.03	   * .data$ether_aromatic),
      b_15 = (-2.23	   * .data$nitrate),
      b_16 = (-2.15	   * .data$nitro),
      b_17 = (-2.14	   * .data$hydroxyl_aromatic),
      b_18 = (-1.03	   * .data$amine_primary),
      b_19 = (-0.849	 * .data$amine_secondary),
      b_20 = (-0.608 	 * .data$amine_tertiary),
      b_21 = (-1.61    * .data$amine_aromatic),
      b_22 = (-4.49    * .data$amide_primary),
      b_23 = (-5.26    * .data$amide_secondary),
      b_24 = (-2.63    * .data$amide_tertiary),
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
        b_32 = (-2.23	  * .data$phosphoric_acids),
        b_33 = (-2.23	  * .data$phosphoric_esters),
        b_34 = (-2.23	  * .data$sulfates),
        b_35 = (-2.23	  * .data$sulfonates),
        b_36 = (-2.23	  * .data$thiols),
        b_37 = (-1.20	  * .data$carbothioesters)
      )
  }
  
  betas %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(log10_P = sum(dplyr::c_across(dplyr::starts_with("b_")))) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-dplyr::starts_with("b_"))
}