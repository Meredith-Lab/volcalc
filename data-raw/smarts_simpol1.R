## code to prepare `smarts` dataset goes here
smarts_simpol1 <- readr::read_csv("data-raw/smarts_simpol1.csv")
usethis::use_data(smarts_simpol1, overwrite = TRUE)
