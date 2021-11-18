test_that("volatility estimate is correct for example compound", {
  ex_fx_df <- get_fx_groups("C16181", "map00361", "data")
  ex_vol_df <- calc_vol(ex_fx_df)
  expect_equal(round(ex_vol_df$log_c, 6), 6.963971)
})
