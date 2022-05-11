test_that("file downloads with correct counts block", {
  save_compound_mol("map00361", path = "tests/testthat/data", compound_id = "C16181", redownload = TRUE)
  ex_file_contents <- readLines("tests/testthat/data/map00361/C16181.mol")
  expect_equal(
    ex_file_contents[4],
    "12 12  0  0  1  0  0  0  0  0999 V2000"
  )
})

test_that("file downloads with formula input with correct counts block", {
  save_compound_mol(pathway_id = "map00361", path = "tests/testthat/data", compound_formula = "C6H4Cl4O2", redownload = TRUE)
  ex_file_contents <- readLines("tests/testthat/data/map00361/C18238.mol")
  expect_equal(
    ex_file_contents[4],
    "12 12  0  0  1  0  0  0  0  0999 V2000"
  )
})

test_that("formula input with no KEGG compound ID match returns error", {
  expect_error(save_compound_mol("map00361", "tests/testthat/data", compound_formula = "C200H200N200"))
})

test_that("formula input with multiple KEGG compound ID matches returns error", {
  expect_error(save_compound_mol("map00623", "tests/testthat/data", compound_formula = "C7H10O5"))
})
