test_that("volatility estimate is correct for example compound for entire workflow", {
  ex_vol_df <- calc_vol("data/C16181.mol")
  expect_equal(round(ex_vol_df$volatility, 6), 6.975349)
})


test_that("returns correct number of columns depending on return arguments", {
  expect_equal(ncol(calc_vol("data/C16181.mol")), 4)
  expect_equal(ncol(calc_vol("data/C16181.mol", return_fx_groups = TRUE)), 43)
  expect_equal(ncol(calc_vol("data/C16181.mol", return_calc_steps = TRUE)), 7)
  expect_equal(ncol(calc_vol("data/C16181.mol", return_fx_groups = TRUE,
                             return_calc_steps = TRUE)), 46)
})

test_that("calc_vol() works with multiple inputs", {
  paths <- c("data/map00361/C00011.mol", "data/map00361/C00042.mol")
  expect_s3_class(calc_vol(paths), "data.frame")
})
  