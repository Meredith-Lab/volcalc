test_that("a functional group count is correct for example compound", {
  ex_df <- get_fx_groups("C16181", "map00361", path = "tests/testthat/data")
  expect_equal(ex_df$hydroxyl_groups, 1)
})
