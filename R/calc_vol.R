#' Calculate volatility estimate for a compound
#'
#' Relative volatility value and category is estimated for specified compound
#' using group contribution methods.
#'
#' @param input a path to a .mol file or a SMILES string.
#' @param from the form of `input`. Either `"mol_path"` or `"smiles"` (default
#'   is `"mol_path"`).
#' @param method the method for calculating estimated volatility. See
#'   [simpol1()] for more details.
#' @param environment the environment for calculating relative volatility
#'   categories. RVI thresholds for low, moderate, and high volatility are as
#'   follows: `"clean"` (clean atmosphere, default) -2, 0, 2; `"polluted"`
#'   (polluted atmosphere) 0, 2, 4; `"soil"` 4, 6, 8.
#' @param return_fx_groups When `TRUE`, includes functional group counts in
#'   final dataframe.
#' @param return_calc_steps When `TRUE`, includes intermediate volatility
#'   calculation steps in final dataframe. See **Details**.
#' 
#' @details \eqn{\textrm{log}_{10}C^\ast} is used for the calculated relative
#'   volatility index (`rvi`). \eqn{\textrm{log}_{10}C^\ast =
#'   \textrm{log}_{10}(PM/RT)} where \eqn{P} is the estimated vapor pressure for
#'   the compound, \eqn{M} is molecular mass of the compound, \eqn{R} is the
#'   universal gas constant, and \eqn{T} is temperature (293.14K or 20ºC).  When
#'   `return_calc_steps = TRUE`, the log of estimated vapor pressure, `log10_P`,
#'   and \eqn{\textrm{log}_{10}(M/RT)}, `log_alpha`, are also returned.
#' 
#'
#' @return a tibble with relative volatility index (`rvi`) and volatility
#'   category (`category`).
#'   
#' @seealso [get_fx_groups()], [simpol1()]
#'
#' @export
#' @examples
#' mol_path <- mol_example("C16181.mol")
#' calc_vol(mol_path)
#' 
#' # Return functional group counts from get_fx_groups()
#' calc_vol(mol_path,  return_fx_groups = TRUE)
#' 
#' # Return intermediate calculations
#' calc_vol(mol_path, return_calc_steps = TRUE)
#' 
calc_vol <-
  function(input, 
           from = c("mol_path", "smiles"),
           method = c("meredith", "simpol1"),
           environment = c("clean", "polluted", "soil"),
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
      #TODO: validate mol files??
      compound_sdf_list <- lapply(input, ChemmineR::read.SDFset)
    }
    
    if(from == "smiles") {
      compound_sdf_list <- lapply(input, ChemmineR::smiles2sdf)
    }
    
    fx_groups_df_list <-
      lapply(compound_sdf_list, get_fx_groups)
    names(fx_groups_df_list) <- input
    fx_groups_df <- 
      #adds column for input named "mol_path" or "smiles"
      dplyr::bind_rows(fx_groups_df_list, .id = {{ from }}) 
    
    # calculate relative volatility & categories from logP
    vol_df <- simpol1(fx_groups_df, meredith = meredith) %>% 
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
      #TODO document log_alpha (need to figure out why it's called that first)
      cols_calc <- c("molecular_weight", "log_alpha", "log10_P")
    }
    
    #return:
    vol_df %>% 
      dplyr::select(dplyr::all_of(c({{ from }}, "formula", "name", "rvi", "category", cols_fx, cols_calc)))
    
  }

