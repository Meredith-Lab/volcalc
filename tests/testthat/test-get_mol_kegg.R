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

test_that("get_mol_kegg errors unless one of compound_id or pathway_id", {
  dir <- withr::local_tempdir()
  expect_error(get_mol_kegg(dir = dir),
               "One of `compound_id` or `pathway_id` are required")
  expect_error(
    get_mol_kegg(compound_ids = "C16181", pathway_ids = "map00361", dir = dir),
    "One of `compound_id` or `pathway_id` are required"
  )
})

test_that("get_mol_kegg checks ID format", {
  dir <- withr::local_tempdir()
  expect_error(
    get_mol_kegg(compound_ids = "hello", dir = dir),
    "Some compound_ids are not in the correct KEGG format"
  )
  expect_error(
    get_mol_kegg(compound_ids = c("banana", "boat"), dir = dir),
    "Some compound_ids are not in the correct KEGG format"
  )
  expect_error(
    get_mol_kegg(compound_ids = c("C16181", "boat"), dir = dir),
    "Some compound_ids are not in the correct KEGG format"
  )
  expect_error(
    get_mol_kegg(pathway_ids = "banana", dir = dir),
    "Some pathway_ids are not in the correct KEGG format"
  )
})

test_that("get_mol_kegg dl correct compound", {
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(compound_ids = "C00083", dir = dir)
  expect_equal(readLines(out$mol_path) %>% head(1), "Malonyl-CoA")
})

test_that("get_mol_kegg works with pathways", {
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(pathway_ids = c("map00253", "map00232"), dir = dir)
  expect_equal(nrow(out), 43)
  expect_equal(unique(out$pathway_id), c("map00253", "map00232"))
})