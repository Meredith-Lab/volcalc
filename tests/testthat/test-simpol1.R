test_that("isoprene is more volatile than glucose", {
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

test_that("changing temperature does something", {
  isoprene <- ChemmineR::read.SDFset(mol_example()[5])
  groups <- get_fx_groups(isoprene)
  c_30 <- simpol1(groups, temp_c = 30)$log10_P
  c_20 <- simpol1(groups, temp_c = 20)$log10_P
  expect_gt(c_30, c_20)
})
