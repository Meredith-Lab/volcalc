test_that("file downloads with correct counts block", {
  save_compound_mol("C16181", "map00361")
  ex_file_contents <- readLines("data/compound_mols/map00361/C16181.mol")
  expect_equal(ex_file_contents[4],
               "12 12  0  0  1  0  0  0  0  0999 V2000")
})
