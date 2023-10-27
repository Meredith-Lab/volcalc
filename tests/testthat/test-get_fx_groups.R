test_that("a functional group count is correct for example compound", {
  sdf <- ChemmineR::read.SDFset("data/C16181.mol")
  ex_df <- get_fx_groups(sdf)
  expect_equal(ex_df$hydroxyl_aliphatic, 1, ignore_attr = TRUE)
})

test_that("error with SDFset with more than one molecule", {
  data(sdfsample,package = "ChemmineR")
  sdfset <- sdfsample
  expect_error(get_fx_groups(sdfset[1:4]),
               "SDFset objects must contain a single molecule only")
})

test_that("SMILES and mol give same results", {
  from_mol <- ChemmineR::read.SDFset("data/C16181.mol")
  from_smiles <- ChemmineR::smiles2sdf("C1(C(C(C(C(C1Cl)Cl)Cl)Cl)Cl)O")
  expect_equal(
    get_fx_groups(from_mol) %>% dplyr::select(-name),
    get_fx_groups(from_smiles) %>% dplyr::select(-name),
    ignore_attr = TRUE
  )
})

test_that("SMARTS strings are correct", {
  # correctness tests compare output of get_fx_groups() to a manually edited
  # .csv file at "tests/testthat/data/test_compounds.csv".  Add compounds to
  # that file with counts of functional groups using the same column names as
  # output by get_fx_groups() to add more tests.
  test_compounds <-
    readr::read_csv("data/test_compounds.csv") %>%
    #for debugging
    # readr::read_csv("tests/testthat/data/test_compounds.csv") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), as.integer))
  test_fx_groups <-
    test_compounds$smiles %>%
    purrr::map(ChemmineR::smiles2sdf) %>%
    purrr::map(get_fx_groups) %>%
    purrr::list_rbind() %>% 
    tibble::add_column(smiles = test_compounds$smiles) %>%
    #remove columns that have only NAs--these don't have SMARTS patterns yet
    dplyr::select(dplyr::where(~!all(is.na(.))))
  
  common_cols <- intersect(
    colnames(test_compounds),
    colnames(test_fx_groups)
  )

  expected <- test_compounds %>% dplyr::select(smiles, dplyr::all_of(common_cols))
  actual   <- test_fx_groups %>% dplyr::select(smiles, dplyr::all_of(common_cols))
  
  # compare but ignore NAs in expected, by just overwriting them with values in actual using rows_patch()
  expect_equal(
    actual,
    dplyr::rows_patch(expected, actual),
    ignore_attr = TRUE
  )
})



