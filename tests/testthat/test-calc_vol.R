test_that("volatility estimate is correct", {
  ex_vol_df <- calc_vol("data/C16181.mol")
  expect_equal(round(ex_vol_df$rvi, 6), 6.975349)
})

test_that("returns correct columns depending on return arguments", {
  just_vol <- calc_vol("data/C16181.mol")
  with_fx <- calc_vol("data/C16181.mol", return_fx_groups = TRUE)
  with_fx_steps <-
    calc_vol("data/C16181.mol",
             return_fx_groups = TRUE,
             return_calc_steps = TRUE)
  expect_setequal(colnames(just_vol),
                  c("mol_path", "formula", "name", "rvi", "category"))
  # just some examples here
  expect_contains(colnames(with_fx),
                  c(colnames(just_vol), "carbons", "carbothioesters", "fluorines"))
  expect_contains(colnames(with_fx_steps),
                  c(colnames(with_fx), "molecular_weight", "log_alpha", "log10_P"))
})

test_that("calc_vol() works with multiple inputs", {
  paths <- c("data/map00361/C00011.mol", "data/map00361/C00042.mol")
  smiles <- c("O=C=O", "C(CC(=O)O)C(=O)O")
  expect_s3_class(calc_vol(paths), "data.frame")
  expect_s3_class(calc_vol(smiles, from = "smiles"), "data.frame")
})

test_that("smiles and .mol give same results", {
  paths <-
    c("data/C16181.mol",
      "data/map00361/C00011.mol",
      "data/map00361/C00042.mol")
  smiles <-
    c("C1(C(C(C(C(C1Cl)Cl)Cl)Cl)Cl)O", "O=C=O", "C(CC(=O)O)C(=O)O")
  expect_equal(
    calc_vol(smiles, from = "smiles") %>% dplyr::select(-name,-smiles),
    calc_vol(paths) %>% dplyr::select(-name,-mol_path)
  )
})

test_that("errors with invalid SMILES", {
  expect_error(calc_vol("hello", from = "smiles"))
})
  
test_that("meredith and original method give different results", {
  #thiol and sulfonate groups, respectively
  # paths <- c(test_path("data/C00409.mol"), test_path("data/C03349.mol"))
  smiles <-
    c("Methanethiol" = "SC",
      "Methyl methanesulfonate" = "COS(=O)(=O)C")
  meredith <- calc_vol(smiles, from = "smiles", method = "meredith")
  simpol   <- calc_vol(smiles, from = "smiles", method = "simpol1")
  expect_true(all(meredith$rvi < simpol$rvi))
})

