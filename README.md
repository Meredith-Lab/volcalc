
<!-- README.md is generated from README.Rmd. Please edit that file -->

# volcalc

<!-- badges: start -->
<!-- badges: end -->

The goal of volcalc is to automate calculating estimates of volatility
for compounds.

## Installation

You can install volcalc from [the GitHub
repository](https://github.com/Meredith-Lab/volcalc) with:

``` r
# install.packages("devtools")
devtools::install_github("Meredith-Lab/volcalc")
```

## Example

This is a basic example which shows you how to get a volatility estimate
for an example compound *beta-2,3,4,5,6-Pentachlorocyclohexanol*. The
KEGG compound identifier for the compound, as found on [the compound’s
KEGG page](https://www.genome.jp/dbget-bin/www_bget?C16181), is
*C16181*. It is part of two molecular pathways; we will use *map00361*
below.

We download the .mol file for the compound from KEGG, count the
occurrences of different functional groups, and use that to estimate
volatility. This calculation uses the SIMPOL approach (Prankow and
Asher, 2008).

``` r
library(volcalc)
save_compound_mol("C16181", "map00361", "data/")
example_compound_fx_groups <- get_fx_groups("C16181", "map00361", "data/")
example_compound_vol <- calc_vol(example_compound_fx_groups)
print(example_compound_vol$log_c)
#> [1] 6.963971
```

This example compound has a volatility around 7. It is in the high
volatility category.
