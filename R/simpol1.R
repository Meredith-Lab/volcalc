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
#' functional groups:
#' 
#' Using the same coefficient as nitrate
#' 
#' - Phosphoric acid
#' - Phosphoric ester
#' - Sulfate
#' - Sulfonate
#' - Thiol
#' 
#' Using the same coefficient as ester
#' - Carbothioester
#' 
#' @note The method described in Pankow & Asher (2008) was developed using data
#'   between 0 °C and 120 °C, so extrapolating beyond that range is not
#'   recommended.
#'
#' @param fx_groups A data.frame or tibble with counts of functional groups
#'   produced by [get_fx_groups()] (or manually, with the same column names).
#' @param meredith Logical; `FALSE`: use the original SIMPOL.1 method. `TRUE`:
#'   use the modified version in Meredith et al. (2023).
#' @param temp_c Numeric; the temperature at which to calculate volatility
#'   estimates in °C.
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
simpol1 <- function(fx_groups, meredith = TRUE, temp_c = 20) {
  #convert C to K
  if (temp_c < 0 | temp_c > 120) {
    warning("Temperatures below 0\u00b0C or above 120\u00b0C extrapolate beyond the range SIMPOL.1 was intended for. \n  Interpret results with caution!")
  }
  temp_k <- temp_c + 273.15
  
  b_k <- calc_bk_simpol1(temp_k)
  
  contributions <- fx_groups %>%
    # assume NAs are 0s for the sake of this calculation
    dplyr::mutate(dplyr::across(dplyr::where(is.integer), function(x)
      ifelse(is.na(x), 0L, x))) %>%
    # TODO: There is also a table of coefs and equation to calculate b_k(T) at
    # different values of T. Why wasn't this implemented?
    
    dplyr::mutate(
      # multiplier for each functional group is volatility contribution to log10P
      # b_k(T)  * v_k,i
      vb_00 = (b_k["b0"]  * 1), #b_0(T) is an intercept/constant.  
      vb_01 = (b_k["b1"]  * .data$carbons),
      vb_02 = (b_k["b2"]  * .data$carbons_asa),
      vb_03 = (b_k["b3"]	 * .data$rings_aromatic),
      vb_04 = (b_k["b4"]  * .data$rings_aliphatic),
      vb_05 = (b_k["b5"]  * .data$carbon_dbl_bonds_aliphatic),
      vb_06 = (b_k["b6"]  * .data$CCCO_aliphatic_ring),
      vb_07 = (b_k["b7"]  * .data$hydroxyl_aliphatic),
      vb_08 = (b_k["b8"]  * .data$aldehydes),
      vb_09 = (b_k["b9"]  * .data$ketones),
      vb_10 = (b_k["b10"]  * .data$carbox_acids),
      vb_11 = (b_k["b11"]  * .data$ester),
      vb_12 = (b_k["b12"]  * .data$ether_alkyl),
      vb_13 = (b_k["b13"]  * .data$ether_alicyclic),
      vb_14 = (b_k["b14"]  * .data$ether_aromatic),
      vb_15 = (b_k["b15"]  * .data$nitrate),
      vb_16 = (b_k["b16"]  * .data$nitro),
      vb_17 = (b_k["b17"]  * .data$hydroxyl_aromatic),
      vb_18 = (b_k["b18"]  * .data$amine_primary),
      vb_19 = (b_k["b19"]	 * .data$amine_secondary),
      vb_20 = (b_k["b20"]	 * .data$amine_tertiary),
      vb_21 = (b_k["b21"]  * .data$amine_aromatic),
      vb_22 = (b_k["b22"]  * .data$amide_primary),
      vb_23 = (b_k["b23"]  * .data$amide_secondary),
      vb_24 = (b_k["b24"]  * .data$amide_tertiary),
      vb_25 = (b_k["b25"]  * .data$carbonylperoxynitrate),
      vb_26 = (b_k["b26"]  * .data$peroxide),
      vb_27 = (b_k["b27"]  * .data$hydroperoxide),
      vb_28 = (b_k["b28"]  * .data$carbonylperoxyacid),
      vb_29 = (b_k["b29"]	 * .data$nitrophenol),
      vb_30 = (b_k["b30"]  * .data$nitroester)
    ) 
  
  if (isTRUE(meredith)) {
    contributions <- contributions %>% 
      dplyr::mutate(
        # # Below are additions from Meredith et al.
        # use the coef for nitrate for the following groups
        vb_32 = (b_k["b15"] * .data$phosphoric_acids),
        vb_33 = (b_k["b15"] * .data$phosphoric_esters),
        vb_34 = (b_k["b15"] * .data$sulfates),
        vb_35 = (b_k["b15"] * .data$sulfonates),
        vb_36 = (b_k["b15"] * .data$thiols),
        # use the coef for ester for carbothioesters
        vb_37 = (b_k["b11"] * .data$carbothioesters)
      )
  }
  
  contributions %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(log10_P = sum(dplyr::c_across(dplyr::starts_with("vb_")))) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(-dplyr::starts_with("vb_"))
}


#from table 5 of Pankow & Asher
calc_bk_simpol1 <- function(temp_k = 293.15) {
  table_5 <- tibble::tribble(
    ~k, ~B1, ~B2, ~B3, ~B4,
    0, -4.26938E+02, 2.89223E-01, 4.42057E-03, 2.92846E-01,
    1, -4.11248E+02, 8.96919E-01, -2.48607E-03, 1.40312E-01,
    2, -1.46442E+02, 1.54528E+00, 1.71021E-03, -2.78291E-01,
    3, 3.50262E+01, -9.20839E-01, 2.24399E-03, -9.36300E-02,
    4, -8.72770E+01, 1.78059E+00, -3.07187E-03, -1.04341E-01,
    5, 5.73335E+00, 1.69764E-02, -6.28957E-04, 7.55434E-03,
    6, -2.61268E+02, -7.63282E-01, -1.68213E-03, 2.89038E-01,
    7, -7.25373E+02, 8.26326E-01, 2.50957E-03, -2.32304E-01,
    8, -7.29501E+02, 9.86017E-01, -2.92664E-03, 1.78077E-01,
    9, -1.37456E+01, 5.23486E-01, 5.50298E-04, -2.76950E-01,
    10, -7.98796E+02, -1.09436E+00, 5.24132E-03, -2.28040E-01,
    11, -3.93345E+02, -9.51778E-01, -2.19071E-03, 3.05843E-01,
    12, -1.44334E+02, -1.85617E+00, -2.37491E-05, 2.88290E-01,
    13, 4.05265E+01, -2.43780E+00, 3.60133E-03, 9.86422E-02,
    14, -7.07406E+01, -1.06674E+00, 3.73104E-03, -1.44003E-01,
    15, -7.83648E+02, -1.03439E+00, -1.07148E-03, 3.15535E-01,
    16, -5.63872E+02, -7.18416E-01, 2.63016E-03, -4.99470E-02,
    17, -4.53961E+02, -3.26105E-01, -1.39780E-04, -3.93916E-02,
    18, 3.71375E+01, -2.66753E+00, 1.01483E-03, 2.14233E-01,
    19, -5.03710E+02, 1.04092E+00, -4.12746E-03, 1.82790E-01,
    20, -3.59763E+01, -4.08458E-01, 1.67264E-03, -9.98919E-02,
    21, -6.09432E+02, 1.50436E+00, -9.09024E-04, -1.35495E-01,
    22, -1.02367E+02, -7.16253E-01, -2.90670E-04, -5.88556E-01,
    23, -1.93802E+03, 6.48262E-01, 1.73245E-03, 3.47940E-02,
    24, -5.26919E+00, 3.06435E-01, 3.25397E-03, -6.81506E-01,
    25, -2.84042E+02, -6.25424E-01, -8.22474E-04, -8.80549E-02,
    26, 1.50093E+02, 2.39875E-02, -3.37969E-03, 1.52789E-02,
    27, -2.03387E+01, -5.48718E+00, 8.39075E-03, 1.07884E-01,
    28, -8.38064E+02, -1.09600E+00, -4.24385E-04, 2.81812E-01,
    29, -5.27934E+01, -4.63689E-01, -5.11647E-03, 3.84965E-01,
    30, -1.61520E+03, 9.01669E-01, 1.44536E-03, 2.66889E-01
  )
  
  b_k_df <- 
    dplyr::mutate(table_5, b_k = .data$B1/temp_k + .data$B2 + .data$B3*temp_k + .data$B4*log(temp_k))
  
  out <- b_k_df$b_k
  names(out) <- paste0("b", b_k_df$k)
  #return:
  out
}