#' A "safe" version of ChemmineR::propOB() that returns NA when OpenBabel
#' produces warnings or errors
#' 
#' Chemminer::propOB() can result in errors from OpenBabel that are difficult to
#' capture because they are not R errors or warnings.  This solution (thanks to
#' Jeron Ooms) runs propOB() in a separate R session and captures stderr in a
#' temp file. By default, if there is anything at all in stderr, the function
#' returns NA, assuming that the results are not reliable.  This behavior can be
#' overriden with `ob_err_pass`, but it is a good idea to inspect warnings and
#' errors with `ob_problems()` and use resulting data at your own risk!
#'
#' @param sdfSet an sdfSet object
#' @param ob_err_pass logical; return results even if there are warnings or errors?
#' 
#' @seealso [ob_problems()]
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun {
#' tmp <- tempdir()
#' sdf <- get_mol_kegg("C00157", dir = tmp)
#' sdf <- ChemmineR::read.SDFset(file.path(tmp, "C00157.mol"))
#' ChemmineR::propOB(sdf)
#' safe_propOB(sdf)
#' }
safe_propOB <- function(sdfSet, ob_err_pass = FALSE) {
  # TODO: make sure these temp dirs get cleaned up.  Use on.exit() or withr
  res <- callr::r(function(sdfSet) ChemmineR::propOB(sdfSet), args = list(sdfSet), stderr = "/tmp/err", stdout = "/tmp/out")
  
  ob_err <- readLines("/tmp/err")
  
  # TODO: consider if you want to distiguish between errors and warnings.  Could be simple as looking for the word "warning" or "error" in ob_err.
  if (length(ob_err) > 0 & isFALSE(ob_err_pass)) {
    warning("ChemmineR::propOB() reported problems, returning NA. View problems with `ob_problems()`.  Disable this behavior with `ob_err_pass = TRUE`")
    res <- NA
    attr(res, "OB_err") <- ob_err
  } 
  res
}

#TODO: document—maybe in same .Rd as safe_propOB
ob_problems <- function(x) {
  attr(x, "OB_err")
}

#TODO: write tests
