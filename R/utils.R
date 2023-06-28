`%+%` <- function(x, y)  mapply(sum, x, y, MoreArgs = list(na.rm = TRUE))

#This only exists to silence a R CMD check note about ChemmineOB being in Imports but no functions being used in volcalc directly.
zzz <- function() {
  ChemmineOB::prop_OB
}
