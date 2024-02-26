#https://ropensci.org/blog/2019/12/08/precompute-vignettes/
knitr::knit("vignettes/kegg-source.Rmd.source", output = "vignettes/kegg.Rmd")
knitr::knit("vignettes/volcalc-source.Rmd.source", output = "vignettes/volcalc.Rmd")