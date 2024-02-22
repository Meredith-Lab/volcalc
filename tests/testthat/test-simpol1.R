test_that("isoprene is more volatile than glucose", {
  skip_if_not_installed("ChemmineOB")
  isoprene <- ChemmineR::read.SDFset(mol_example()[5])
  glucose  <- ChemmineR::read.SDFset(mol_example()[1])
  
  log10P_together <- simpol1(dplyr::bind_rows(
    get_fx_groups(isoprene),
    get_fx_groups(glucose)
  ))$log10_P
  
  log10P_separate <- c(
    simpol1(get_fx_groups(isoprene))$log10_P,
    simpol1(get_fx_groups(glucose))$log10_P
  )
  
  expect_gt(log10P_together[1], log10P_together[2])
  expect_equal(log10P_together, log10P_separate)
})
