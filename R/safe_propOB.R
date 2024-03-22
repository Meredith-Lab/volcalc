#' "Safe" version of ChemmineR::propOB()
#' 
#' A "safe" version of [ChemmineR::propOB()] that returns `NA`s when OpenBabel
#' produces warnings or errors.
#'
#' [ChemmineR::propOB()] can result in errors from OpenBabel that are difficult
#' to capture because they are not R errors or warnings (i.e. they are not
#' captured by [capture.output()] or [sink()]). This solution (thanks to Jeron
#' Ooms) runs [ChemmineR::propOB()] in a separate R session and captures
#' `stderr` in a temp file. By default, if there is anything at all in `stderr`,
#' all values returned are `NA`s, assuming that the results are not reliable.
#' This behavior can be overridden with `ob_err_pass`, but it is a good idea to
#' inspect warnings and errors with [ob_problems()] and use resulting data at
#' your own risk!
#'
#' @param sdfSet an `sdfSet` object
#' @param ob_err_pass logical; return results even if there are warnings or
#'   errors? Defaults to `FALSE`
#'
#' @seealso [ob_problems()]
#'
#' @return A data frame with any messages from OpenBabel stored as the
#'   attribute `"OB_stderr"`.  See [ChemmineR::propOB()] for additional details.
#' @export
#'
#' @examples
#' \dontrun{
#' tmp <- tempdir()
#' sdf <- get_mol_kegg("C00157", dir = tmp)
#' sdf <- ChemmineR::read.SDFset(file.path(tmp, "C00157.mol"))
#' prop <- ChemmineR::propOB(sdf)
#' prop_safe <- safe_propOB(sdf)
#' ob_problems(prop_safe)
#' }
safe_propOB <- function(sdfSet, ob_err_pass = FALSE) {
  #TODO: consider two separate args:
   # - one to capture the errors/warnings but ignore them and report the values anyway.
   # - one passthrough to propOB()
  # worth it? Maybe not.  Probably should do benchmarking
  temp_err <- tempfile()
  on.exit(file.remove(temp_err))
  res <- 
    callr::r(function(sdfSet) ChemmineR::propOB(sdfSet), args = list(sdfSet), stderr = temp_err)
  ob_err <- readLines(temp_err)
  
  # TODO: consider if you want to distinguish between errors and warnings.  Could
  # be simple as looking for the word "warning" or "error" in ob_err. Perhaps I
  # should catalog the different kinds of warnings/errors we get?
  
  if (length(ob_err) > 0 & isFALSE(ob_err_pass)) {
    #TODO think about how these warnings will be displayed when safe_propOB is run on 100+ SDFs.  Should they be combined?  One warning for every problem? Should the warning mention the compound name or SMILES or something?
    cli::cli_warn(c(
      "!" ="{.fn ChemmineR::propOB} reported problems, returning {.code NA}.",
      "i" = "View problems with {.fn ob_problems}.",
      "i" = "Disable this behavior with {.code ob_err_pass = TRUE}."
    ))
    res[ , colnames(res) != "title"] <- NA
  }
  attr(res, "OB_stderr") <- ob_err
  res
}

#' @param x data frame created by [safe_propOB()]
#' @rdname safe_propOB
#' @export
ob_problems <- function(x) {
  OB_stderr <- attr(x, "OB_stderr")
  cat(OB_stderr, sep = "\n")
  return(invisible(OB_stderr))
}

#TODO: write tests

## Some benchmarking
# sdfs <- fs::dir_ls("tests/testthat/data/", glob = "*.mol", recurse = TRUE) %>%
#   purrr::map(ChemmineR::read.SDFset)
# 
# library(tictoc)
# library(purrr)
# # capture errors and return NAs
# tic()
# x <- purrr::map(sdfs, safe_propOB)
# toc() # 242.718 sec elapsed
# list_rbind(x)
# 
# # vs. regular propOB()
# tic()
# x2 <- purrr::map(sdfs, ChemmineR::propOB)
# toc() # 0.701 sec elapsed



