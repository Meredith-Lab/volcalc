test_that("volatility estimate is correct for example compound for entire workflow", {
  ex_vol_df <- calc_vol("data/C16181.mol")
  expect_equal(round(ex_vol_df$rvi, 6), 6.975349)
})


test_that("returns correct number of columns depending on return arguments", {
  #TODO re-write test to check for number of columns *more* than baseline so adding functional groups doesn't break tests
  expect_equal(ncol(calc_vol("data/C16181.mol")), 5)
  expect_equal(ncol(calc_vol("data/C16181.mol", return_fx_groups = TRUE)), 51)
  expect_equal(ncol(calc_vol("data/C16181.mol", return_calc_steps = TRUE)), 8)
  expect_equal(ncol(calc_vol("data/C16181.mol", return_fx_groups = TRUE,
                             return_calc_steps = TRUE)), 54)
})

test_that("calc_vol() works with multiple inputs", {
  paths <- c("data/map00361/C00011.mol", "data/map00361/C00042.mol")
  smiles <- c("O=C=O", "C(CC(=O)O)C(=O)O")
  expect_s3_class(calc_vol(paths), "data.frame")
  expect_s3_class(calc_vol(smiles, from = "smiles"), "data.frame")
})

test_that("smiles and .mol give same results", {
  paths <- c("data/C16181.mol", "data/map00361/C00011.mol", "data/map00361/C00042.mol")
  smiles <- c("C1(C(C(C(C(C1Cl)Cl)Cl)Cl)Cl)O", "O=C=O", "C(CC(=O)O)C(=O)O")
  expect_equal(
    calc_vol(smiles, from = "smiles") %>% dplyr::select(-name, -smiles),
    calc_vol(paths) %>% dplyr::select(-name, -mol_path)
  )
})

test_that("errors with invalid SMILES", {
  expect_error(calc_vol("hello", from = "smiles"))
})
