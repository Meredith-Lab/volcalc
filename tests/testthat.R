library(testthat)
library(volcalc)

testthat::skip_on_os("windows"){
  ## Skip on Windows because ChemmineOB package not available; 
  ## tests fail
  test_check("volcalc")  
}

