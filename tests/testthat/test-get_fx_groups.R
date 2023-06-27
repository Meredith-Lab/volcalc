test_that("a functional group count is correct for example compound", {
  sdf <- ChemmineR::read.SDFset("data/C16181.mol")
  ex_df <- get_fx_groups(sdf)
  expect_equal(ex_df$hydroxyl_groups, 1)
})

#' TODO:
#' - Additional correctness tests
#' - Test that compound name is read correctly or is NA when empty
