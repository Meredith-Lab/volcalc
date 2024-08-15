#' Calculate volatility estimate for a compound
#'
#' Estimate relative volatility index value and category for specified compounds
#' using group contribution methods.
#'
#' @param input A path to a .mol file or a SMILES string.
#' @param from The form of `input`. Either `"mol_path"` (default) or `"smiles"`.
#' @param method The method for calculating estimated volatility. See
#'   [simpol1()] for more details.
#' @param temp_c numeric; For methods that allow it, specify temperature (in
#'   degrees C) at which log10(P) estimates are calculated.
#' @param environment The environment for calculating relative volatility
#'   categories. RVI thresholds for low, moderate, and high volatility are as
#'   follows: `"clean"` (clean atmosphere, default) -2, 0, 2; `"polluted"`
#'   (polluted atmosphere) 0, 2, 4; `"soil"` 4, 6, 8. For more information about
#'   these thresholds see Meredith et al. (2023) and Donahue et al. (2006).
#' @param validate logical; if `TRUE` (default), results are checked for
#'   possible errors in parsing by Open Babel and `NA`s are returned if possible
#'   errors are found. Setting to `FALSE` bypasses these checks—use at your own
#'   risk! Validation is not available on Windows. See **Details** of
#'   [get_fx_groups()] for more information.
#' @param return_fx_groups When `TRUE`, the returned tibble includes functional
#'   group counts.
#' @param return_calc_steps When `TRUE`, the returned tibble includes
#'   intermediate volatility calculations. See **Details**.
#'
#' @details \eqn{\textrm{log}_{10}C^\ast} is used for the calculated relative
#'   volatility index (`rvi`). \eqn{\textrm{log}_{10}C^\ast =
#'   \textrm{log}_{10}(PM/RT)} where \eqn{P} is the estimated vapor pressure for
#'   the compound, \eqn{M} is molecular mass of the compound, \eqn{R} is the
#'   universal gas constant, and \eqn{T} is temperature (293.14K or 20ºC).  When
#'   `return_calc_steps = TRUE`, the log of estimated vapor pressure, `log10_P`,
#'   and \eqn{\textrm{log}_{10}(M/RT)}, `log_alpha`, are also returned.
#'
#' @references Donahue, N.M., Robinson, A.L., Stanier, C.O., Pandis, S.N., 2006.
#' Coupled Partitioning, Dilution, and Chemical Aging of Semivolatile Organics.
#' Environ. Sci. Technol. 40, 2635–2643. \doi{10.1021/es052297c}
#'
#' Meredith L, Ledford S, Riemer K, Geffre P, Graves K, Honeker L, LeBauer D,
#' Tfaily M, Krechmer J. 2023. Automating methods for estimating metabolite
#' volatility. Frontiers in Microbiology. \doi{10.3389/fmicb.2023.1267234}
#'
#' @returns A tibble with relative volatility index (`rvi`) and volatility
#'   category (`category`).
#'
#' @seealso [get_fx_groups()], [simpol1()]
#'
#' @export
#' @examples
#' \donttest{
#' mol_paths <- mol_example()
#' calc_vol(mol_paths)
#'
#' # Return functional group counts from get_fx_groups()
#' calc_vol(mol_paths,  return_fx_groups = TRUE)
#'
#' # Return intermediate calculations
#' calc_vol(mol_paths, return_calc_steps = TRUE)
#' }
calc_vol <-
  function(input, 
           from = c("mol_path", "smiles"),
           method = c("meredith", "simpol1"),
           temp_c = 20,
           environment = c("clean", "polluted", "soil"),
           validate = TRUE,
           return_fx_groups = FALSE,
           return_calc_steps = FALSE) {
    
    from <- match.arg(from)
    
    # logic here will likely need to change if new method functions are added
    method <- match.arg(method)
    
    if (method == "meredith") {
      meredith <- TRUE
    } else {
      meredith <- FALSE
    }
    environment <- match.arg(environment)
    
    cutoffs <- switch(
      environment,
      "clean" = c(-Inf, -2, 0, 2, Inf),
      "polluted" = c(-Inf, 0, 2, 4, Inf),
      "soil" = c(-Inf, 4, 6, 8, Inf)
    )

    if(from == "mol_path") {
      compound_sdf_list <- lapply(input, ChemmineR::read.SDFset)
    }
    
    if(from == "smiles") {
      # a way of converting into a list of *named* vectors, so smiles2sdf() runs on named vectors as input
      input <- split(input, seq_along(input))
      compound_sdf_list <- lapply(input, ChemmineR::smiles2sdf)
    }
    fx_groups_df_list <-
      lapply(compound_sdf_list, get_fx_groups, validate = validate)
    names(fx_groups_df_list) <- input
    fx_groups_df <- 
      #adds column for input named "mol_path" or "smiles"
      dplyr::bind_rows(fx_groups_df_list, .id = {{ from }}) 
    
    # calculate relative volatility & categories from logP
    vol_df <- simpol1(fx_groups_df, meredith = meredith, temp_c = temp_c) %>% 
      dplyr::mutate(
        # mass is converted from grams to micrograms
        # 0.0000821 is universal gas constant
        # 293.15 is temperature in Kelvins (20ºC)
        log_alpha = log10((1000000 * .data$molecular_weight) / (0.0000821 * 293.15)),
        rvi = .data$log_alpha + .data$log10_P, 
        category = cut(
          .data$rvi,
          breaks = cutoffs,
          labels = c("non-volatile", "low", "moderate", "high"),
          right = FALSE
        )
      )
    
    # wrangle output
    cols_fx <- NULL
    cols_calc <- NULL
    if (isTRUE(return_fx_groups)) {
      cols_fx <- 
        colnames(fx_groups_df)[!colnames(fx_groups_df) %in% c("formula", "name", "molecular_weight")]
    }
    if (isTRUE(return_calc_steps)) {
      cols_calc <- c("molecular_weight", "log_alpha", "log10_P")
    }
    
    #return:
    vol_df %>% 
      dplyr::select(dplyr::all_of(c(
        {{ from }}, "formula", "name", "rvi", "category", cols_fx, cols_calc)
      ))
    
  }

