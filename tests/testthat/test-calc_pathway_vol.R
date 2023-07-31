test_that("average volatility is correct for all compounds in a pathway ", {
  skip_on_cran()
  ex_df <- calc_pathway_vol("map00361", path = withr::local_tempdir())
  expect_equal(round(mean(ex_df$volatility), 3), 4.597)
})
