
<!-- README.md is generated from README.Rmd. Please edit that file -->

# volcalc

<!-- badges: start -->
<!-- badges: end -->

The goal of volcalc is to automate calculating estimates of volatility
for compounds.

## Installing and loading package

You can install volcalc from [the GitHub
repository](https://github.com/Meredith-Lab/volcalc) with:

``` r
# install.packages("devtools")
devtools::install_github("Meredith-Lab/volcalc")
```

And then use the package with:

``` r
library(volcalc)
```

## Single compound usage

This is a basic example which shows you how to get a volatility estimate
for an example compound *beta-2,3,4,5,6-Pentachlorocyclohexanol*. The
KEGG compound identifier for the compound, as found on [the compound’s
KEGG page](https://www.genome.jp/dbget-bin/www_bget?C16181), is
*C16181*. It is part of two molecular pathways; we will use *map00361*
below.

#### Single function approach

``` r
calc_vol("map00361", "data", compound_id = "C16181")
#>       pathway compound  formula                                   name    log_c
#> CMP1 map00361   C16181 C6H7Cl5O beta-2,3,4,5,6-Pentachlorocyclohexanol 6.963971
#>      volatility
#> CMP1       high
```

This returns a dataframe with columns specifying general info about the
compound, and the compound’s calculated volatility and corresponding
volatility category. The functional group counts underlying the
volatility can be additionally returned with `return_fx_groups = TRUE`,
and the intermediate calculation steps with `return_calc_steps = TRUE`.
A list of all possible dataframe columns is included below.

The compound can alternatively be specified with its chemical formula
using the `compound_formula` argument instead of `compound_id` as in the
example.

#### Multiple function approach

This breaks the steps done by `calc_vol` into three parts: 1) download
the compound’s .mol file from KEGG, 2) count occurrences of different
functional groups, and 3) estimate volatility. This calculation uses the
SIMPOL approach (Prankow and Asher, 2008).

``` r
save_compound_mol("map00361", "data/", compound_id = "C16181")
example_compound_fx_groups <- get_fx_groups("C16181", "map00361", "data/")
example_compound_vol <- calc_vol("map00361", "data/", compound_id = "C16181", fx_groups_df = example_compound_fx_groups)
print(example_compound_vol$log_c)
#> [1] 6.963971
```

This example compound has a volatility around 7. It is in the high
volatility category.

## Multiple compounds from a pathway usage

A dataframe with volatility estimates for all compounds in a chosen
pathway can be returned as below.

``` r
example_pathway_vol <- calc_pathway_vol("map00361", "data/")
print(example_pathway_vol[1,])
#>       pathway compound formula name    log_c volatility
#> CMP1 map00361   C00011     CO2 CO2; 7.912336       high
```

## Dataframe columns

Basic compound information

-   pathway: KEGG pathway identifier
-   compound: KEGG compound identifier
-   formula: compound chemical formula
-   name: compound name
-   mass: compound mass

Counted functional groups and atoms

-   carbons  
-   ketones  
-   aldehydes
-   hydroxyl\_groups
-   carbox\_acids  
-   peroxide
-   hydroperoxide  
-   nitrate  
-   nitro  
-   carbon\_dbl\_bonds
-   rings  
-   rings\_aromatic  
-   phenol  
-   nitrophenol  
-   nitroester  
-   ester  
-   ether\_alicyclic
-   ether\_aromatic  
-   amine\_primary  
-   amine\_secondary
-   amine\_tertiary  
-   amine\_aromatic  
-   amines  
-   amides  
-   phosphoric\_acid
-   phosphoric\_ester
-   sulfate  
-   sulfonate
-   thiol  
-   carbothioester  
-   oxygens  
-   chlorines
-   nitrogens
-   sulfurs  
-   phosphoruses  
-   bromines
-   iodines  
-   fluorines

Volatility calculation steps

-   log\_alpha: intermediate step
-   log\_Sum: intermediate step
-   log\_c: estimated volatility
-   volatility: volatility category, where values less than 0 are
    “none”, values between 0 and 2 are “moderate”, and values above 2
    are “high”
