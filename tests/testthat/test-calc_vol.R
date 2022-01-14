test_that("volatility estimate is correct for example compound for entire workflow", {
  ex_vol_df <- calc_vol("map00361", "data", compound_id = "C16181")
  expect_equal(round(ex_vol_df$log_c, 6), 6.963971)
})

test_that("volatility estimate is correct for example compound from functional groups dataframe", {
  ex_fx_df <- get_fx_groups("C16181", "map00361", "data")
  ex_vol_df <- calc_vol("map00361", "data", compound_id = "C16181", save_file = FALSE,
                        get_groups = FALSE, fx_groups_df = ex_fx_df)
  expect_equal(round(ex_vol_df$log_c, 6), 6.963971)
})

test_that("returns error if no functional groups dataframe and no saving file and getting functional groups dataframe", {
  expect_error(calc_vol("map00361", "data", compound_id = "C16181", save_file = FALSE,
                        get_groups = FALSE, fx_groups_df = NULL))
})
