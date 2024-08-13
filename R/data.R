#' Search patterns used for SIMPOL.1 functional groups
#'
#' This dataframe documents how functional groups for the SIMPOL.1 and Meredith
#' et al. method are defined using SMARTS strings or `ChemmineR` functions.
#'
#' @format
#' \describe{
#'   \item{method}{Either "simpol1" for functional groups only used with the SIMPOL.1 method, or "meredith" for additional groups used in the Meredith et al. method.}
#'   \item{functional_groups}{These correspond to matching column names in the results of [get_fx_groups()].}
#'   \item{description}{Functional group description from Table 5 of Pankow & Asher (2008)}
#'   \item{smarts}{SMARTS strings used to capture groups, when applicable}
#'   \item{fun}{The function used to capture the functional group.  When `smarts` is not `NA`, this is always "[ChemmineR::smartsSearchOB]".  Other groups are captured with other `ChemmineR` functions or as calculations using other functional groups.}
#'   \item{notes}{Notes including how any functional group counts are corrected when there is overlap.  E.g. when one SMARTS pattern is a subset of another pattern, but the two groups are counted separately without overlap in the SIMPOL.1 method.}
#' }
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
"smarts_simpol1"