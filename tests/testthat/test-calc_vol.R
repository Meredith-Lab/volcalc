test_that("volatility estimate is correct for example compound for entire workflow", {
  ex_vol_df <- calc_vol("map00361", path = "tests/testthat/data", compound_id = "C16181")
  expect_equal(round(ex_vol_df$volatility, 6), 6.963571)
})

test_that("volatility estimate is correct for example compound from functional groups dataframe", {
  ex_fx_df <- get_fx_groups("C16181", "map00361", path = "tests/testthat/data")
  ex_vol_df <- calc_vol("map00361", path = "tests/testthat/data", compound_id = "C16181", save_file = FALSE,
                        get_groups = FALSE, fx_groups_df = ex_fx_df)
  expect_equal(round(ex_vol_df$volatility, 6), 6.963571)
})

test_that("returns error if no functional groups dataframe and no saving file and getting functional groups dataframe", {
  expect_error(calc_vol("map00361", path = "tests/testthat/data", compound_id = "C16181", save_file = FALSE,
                        get_groups = FALSE, fx_groups_df = NULL))
})

test_that("returns correct number of columns depending on return arguments", {
  expect_equal(ncol(calc_vol("map00361", path = "tests/testthat/data", compound_id = "C16181")), 6)
  expect_equal(ncol(calc_vol("map00361", path = "tests/testthat/data", compound_id = "C16181", return_fx_groups = TRUE)), 44)
  expect_equal(ncol(calc_vol("map00361", path = "tests/testthat/data", compound_id = "C16181", return_calc_steps = TRUE)), 9)
  expect_equal(ncol(calc_vol("map00361", path = "tests/testthat/data", compound_id = "C16181", return_fx_groups = TRUE, return_calc_steps = TRUE)), 47)
})
