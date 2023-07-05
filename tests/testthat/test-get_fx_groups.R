test_that("a functional group count is correct for example compound", {
  sdf <- ChemmineR::read.SDFset("data/C16181.mol")
  ex_df <- get_fx_groups(sdf)
  expect_equal(ex_df$hydroxyl_groups, 1)
})

test_that("error with SDFset with more than one molecule", {
  data(sdfsample,package = "ChemmineR")
  sdfset <- sdfsample
  expect_error(get_fx_groups(sdfset[1:4]))
})
#' TODO:
#' - Additional correctness tests
#' - Test that compound name is read correctly or is NA when empty

test_that("get_fx_groups() distinguishes between ROOR, ROOH, and ROH", {
  skip("no pattern for hydroperoxide yet")
  #https://en.wikipedia.org/wiki/Tert-Butyl_hydroperoxide
  hydroperoxide <- ChemmineR::smiles2sdf("CC(C)(C)OO")
  
  #https://en.wikipedia.org/wiki/Dicumyl_peroxide
  peroxide <- ChemmineR::smiles2sdf("CC(C)(C1=CC=CC=C1)OOC(C)(C)C2=CC=CC=C2") 
  
  #https://en.wikipedia.org/wiki/Glucose
  alcohol <- ChemmineR::smiles2sdf("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O")
  
  rooh <- get_fx_groups(hydroperoxide)
  roor <- get_fx_groups(peroxide)
  roh  <- get_fx_groups(alcohol)
  
  expect_equal(rooh$hydroperoxide, 1)
  expect_equal(rooh$peroxide, 0) #don't double-count peroxides and hydroperoxides
  expect_equal(rooh$hydroxyl_groups, 0)
  
  expect_equal(roor$peroxide, 1)
  expect_equal(roor$hydroxyl_groups, 0)
  expect_equal(roor$hydroperoxide, 0)
  
  expect_equal(roh$hydroxyl_groups, 5)
  expect_equal(roh$peroxide, 0)
  expect_equal(roh$hydroperoxide, 0)
  
})
