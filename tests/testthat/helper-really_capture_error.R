# There is some debate on where to put these kinds of functions
# (https://blog.r-hub.io/2020/11/18/testthat-utility-belt/).  I chose here
# because I don't want to deal with processx being in R/ but not in Imports

# This was modified from code written by Gábor Csárdi.  The purpose is to
# capture OpenBabel system library errors that are not captured by
# capture.output() or sink() despite them being printed to the R console.  This
# is a "very nasty trick" according to Gábor so take care in using this
really_capture_error <- function(expr) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  orig <- processx::conn_set_stderr(
    processx::conn_create_file(tmp, write = TRUE),
    drop = FALSE
  )
  on.exit(processx::conn_set_stderr(orig), add = TRUE)
  expr
  stderr <- readLines(tmp)
  if(length(stderr) != 0) {
    stop(stderr)
  }
}

