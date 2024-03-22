test_that("get_mol_kegg returns tibble", {
  skip_if_offline()
  skip_on_cran()
  
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(compound_ids = "C16181", dir = dir)
  expect_s3_class(out, "tbl_df")
  expect_equal(out$mol_path, fs::path(dir, "C16181", ext = "mol"))
})

test_that("get_mol_kegg writes files", {
  skip_if_offline()
  skip_on_cran()
  
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(compound_ids = "C16181", dir = dir)
  expect_true(fs::file_exists(fs::path(out$mol_path)))
})

test_that("get_mol_kegg errors unless one of compound_ids or pathway_ids", {
  skip_if_offline()
  skip_on_cran()
  
  dir <- withr::local_tempdir()
  expect_error(get_mol_kegg(dir = dir),
               "One of `compound_ids` or `pathway_ids` are required")
  expect_error(
    get_mol_kegg(compound_ids = "C16181", pathway_ids = "map00361", dir = dir),
    "One of `compound_ids` or `pathway_ids` are required"
  )
})

test_that("get_mol_kegg checks ID format", {
  skip_if_offline()
  skip_on_cran()
  
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
  skip_if_offline()
  skip_on_cran()
  
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(compound_ids = "C00083", dir = dir)
  expect_equal(readLines(out$mol_path) %>% head(1), "Malonyl-CoA")
})

test_that("get_mol_kegg works with pathways", {
  skip_if_offline()
  skip_on_cran()
    
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(pathway_ids = c("map00253", "map00232"), dir = dir)
  expect_equal(nrow(out), 43)
  expect_equal(unique(out$pathway_id), c("map00253", "map00232"))
})

test_that(".mol files are correctly formed", {
  skip_if_offline()
  skip_on_cran()
  skip_on_os("windows") # really_capture_error() errors on windows
  
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(compound_ids = "C00083", dir = dir)
  sdf <- ChemmineR::read.SDFset(out$mol_path)
  # sdf <- ChemmineR::read.SDFset("tests/testthat/data/map00361/C00042.mol")
  expect_no_error(
    really_capture_error(x <- ChemmineR::propOB(sdf))
  )
  
})

test_that("file downloads with correct counts block", {
  skip_on_cran()
  skip_if_offline()
  
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(compound_ids = "C16181", dir = dir)
  ex_file_contents <- readLines(out$mol_path)
  expect_equal(
    ex_file_contents[4],
    "12 12  0  0  1  0  0  0  0  0999 V2000"
  )
})

test_that("works with pathway modules", {
  skip_on_cran()
  skip_if_offline()
  
  dir <- withr::local_tempdir()
  out <- get_mol_kegg(pathway_ids = "M00082", dir = dir)
  expect_equal(nrow(out), 5)
  expect_true(all(file.exists(out$mol_path)))
})

