## code to prepare `smarts` dataset goes here
smarts_simpol1 <- readr::read_csv("data-raw/smarts_simpol1.csv")

#create user-facing data.frame
usethis::use_data(smarts_simpol1, overwrite = TRUE)

#create internal named list with just SMARTS strings
just_smarts_simpol1 <- 
  smarts_simpol1 %>% 
  dplyr::filter(!is.na(smarts)) 
smarts_patterns_simpol1 <- as.list(just_smarts_simpol1$smarts)
names(smarts_patterns_simpol1) <- just_smarts_simpol1$functional_group

usethis::use_data(smarts_patterns_simpol1, internal = TRUE, overwrite = TRUE)
