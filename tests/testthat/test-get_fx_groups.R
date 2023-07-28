test_that("a functional group count is correct for example compound", {
  ex_df <- get_fx_groups(compound_id = "C16181", path = "data")
  expect_equal(ex_df$hydroxyl_groups, 1)
})

test_that("functional group count is correct for example compound with specified pathway", {
  ex_df <- get_fx_groups(compound_id = "C16181", pathway_id = "map00361",
                         path = "data")
  expect_equal(ex_df$hydroxyl_groups, 1)
})

test_that("correct number of rings", {
  ex_df <- get_fx_groups(compound_id = "C07481", path = "data")
  bpa <- get_fx_groups(compound_id = "C13624", path = "data")
  expect_equal(ex_df$rings_aromatic, 2)
  expect_equal(ex_df$rings, 0)
  expect_equal(bpa$rings_aromatic, 2)
  expect_equal(bpa$rings, 0)
})

test_that("correct number of aromatic hydroxyl", {
  bpa <- get_fx_groups(compound_id = "C13624", path = "data")
  expect_equal(bpa$hydroxyl_groups, 0)
  expect_equal(bpa$hydroxyl_aromatic, 2)
  ldopa <- get_fx_groups(compound_id = "C00355", path = "data")
  expect_equal(ldopa$hydroxyl_groups, 0)
  expect_equal(ldopa$hydroxyl_aromatic, 2)
})