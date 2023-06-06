test_that("a functional group count is correct for example compound", {
  ex_df <- get_fx_groups(compound_id = "C16181", path = "data")
  expect_equal(ex_df$hydroxyl_groups, 1)
})

test_that("functional group count is correct for example compound with specified pathway", {
  ex_df <- get_fx_groups(compound_id = "C16181", pathway_id = "map00361",
                         path = "data")
  expect_equal(ex_df$hydroxyl_groups, 1)
})
