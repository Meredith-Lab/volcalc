test_that("average volatility is correct for all compounds in a pathway ", {
  ex_df <- calc_pathway_vol("map00361", "data/")
  expect_equal(round(mean(ex_df$volatility), 3), 4.608)
})
