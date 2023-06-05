# These test that the download from the KEGG API is working.  Skipped on CRAN to
# avoid problems if API is down.
test_that("file downloads with correct counts block", {
  skip_on_cran()
  
  out <- withr::local_tempdir()
  save_compound_mol(compound_id = "C16181", path = out)
  ex_file_contents <- readLines(file.path(out, "C16181.mol"))
  expect_equal(
    ex_file_contents[4],
    "12 12  0  0  1  0  0  0  0  0999 V2000"
  )
})

test_that("file downloads into pathway subfolder with correct counts block", {
  skip_on_cran()
  
  out <- withr::local_tempdir()
  save_compound_mol(compound_id = "C16181",
                    pathway_id = "map00361",
                    path = out)
  ex_file_contents <- readLines(file.path(out, "map00361/C16181.mol"))
  expect_equal(
    ex_file_contents[4],
    "12 12  0  0  1  0  0  0  0  0999 V2000"
  )
})

test_that("file downloads with formula input with correct counts block", {
  skip_on_cran()
  
  out <- withr::local_tempdir()
  save_compound_mol(compound_formula = "C6H4Cl4O2",
                    pathway_id = "map00361",
                    path = out)
  ex_file_contents <- readLines(file.path(out, "map00361/C18238.mol"))
  expect_equal(
    ex_file_contents[4],
    "12 12  0  0  1  0  0  0  0  0999 V2000"
  )
})

test_that("formula input with no KEGG compound ID match returns error", {
  skip_on_cran()
  
  expect_error(
    save_compound_mol(
      compound_formula = "C200H200N200",
      pathway_id = "map00361",
      path = withr::local_tempdir()
    )
  )
})

test_that("formula input with multiple KEGG compound ID matches returns error", {
  skip_on_cran()
  
  expect_error(
    save_compound_mol(
      compound_formula = "C7H10O5",
      pathway_id = "map00623",
      path = withr::local_tempdir()
    )
  )
})
