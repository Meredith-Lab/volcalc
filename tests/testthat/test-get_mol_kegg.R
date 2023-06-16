test_that("get_mol_kegg returns tibble", {
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(compound_ids = "C16181", dir = dir)
  expect_s3_class(out, "tbl_df")
  expect_equal(out$mol_path, fs::path(dir, "C16181", ext = "mol"))
})

test_that("get_mol_kegg writes files", {
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(compound_ids = "C16181", dir = dir)
  expect_true(fs::file_exists(fs::path(out$mol_path)))
})