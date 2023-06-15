#This only exists to silence a R CMD check note about ChemmineOB being in Imports but no functions being used in volcalc directly.
dummy <- function() {
  ChemmineOB::prop_OB
}