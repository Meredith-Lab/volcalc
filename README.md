
<!-- README.md is generated from README.Rmd. Please edit that file -->

# volcalc <a href="https://meredith-lab.github.io/volcalc/"><img src="man/figures/logo.PNG" alt="volcalc website" align="right" height="120"/></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/Meredith-Lab/volcalc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Meredith-Lab/volcalc/actions/workflows/R-CMD-check.yaml)
[![latest-DOI](https://zenodo.org/badge/425022983.svg)](https://zenodo.org/badge/latestdoi/425022983)
[![manuscript-DOI](https://img.shields.io/badge/DOI-10.3389/fmicb.2023.1267234-32a859.svg)](https://doi.org/10.3389/fmicb.2023.1267234)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Codecov test
coverage](https://codecov.io/gh/Meredith-Lab/volcalc/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Meredith-Lab/volcalc?branch=master)
[![volcalc status
badge](https://cct-datascience.r-universe.dev/badges/volcalc)](https://cct-datascience.r-universe.dev/volcalc)

<!-- badges: end -->

## Overview

The `volcalc` package allows you to automate calculating estimates of
volatility for chemical compounds.

> [!WARNING]
> `volcalc` is a work in progress---use at your own risk!

For a bit of a road map of where development is headed, see our
[proposal](https://cct-datascience.github.io/volcalc-isc-proposal/) for
the R Consortium grant.

`volcalc` is designed to support “group contribution” methods for
estimating volatility that rely on molecular properties such as
molecular weight, numbers of certain atoms, and counts of certain
functional groups. Currently, the only methods implemented are SIMPOL.1
(Pankow & Asher 2008) and a modified version used in Meredith et
al. (2023).

`volcalc` works with either .mol files or
[SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system)
strings as input, and supports downloading .mol files directly from
[KEGG](https://www.kegg.jp/).

## Installation

You can install the development version of `volcalc` from GitHub with

``` r
# install.packages("pak")
pak::pkg_install("Meredith-Lab/volcalc")
```

Or from r-universe with

``` r
install.packages("volcalc", repos = c("https://cct-datascience.r-universe.dev", getOption("repos")))
```

You can install the ‘legacy’ version used in Meredith et al. (2023) with

``` r
pak::pkg_install("Meredith-Lab/volcalc@v1.0.2")
```

Installation of `volcalc` requires the system libraries
[OpenBabel](https://openbabel.org/wiki/Main_Page) and Eigen3
(requirements of the `ChemmineOB` package, which `volcalc` depends on).
`pak` will take care of the installation of these libraries for you on
some systems, but you may need to install them manually.

For macOS, they can be installed via homebrew by running the following
shell command:

``` bash
brew install open-babel
```

For Ubuntu Linux:

``` bash
sudo apt-get install libopenbabel-dev
sudo apt-get install libeigen3-dev
```

For windows, `OpenBabel` is included in the `ChemmineOB` binary and does
not need to be installed separately.

For other installation options see the [OpenBabel
documentation](https://openbabel.org/docs/dev/Installation/install.html)
and `ChemmineOB` [install
guide](https://github.com/girke-lab/ChemmineOB/blob/master/INSTALL)

## Basic Usage

This is a basic example which shows you how to get an estimated relative
volatility index (`rvi`) for two example compounds
*beta-2,3,4,5,6-Pentachlorocyclohexanol*, and *Succinate*. The KEGG
compound identifiers for the compounds, as found on [the compound’s KEGG
page](https://www.genome.jp/dbget-bin/www_bget?C16181), are C16181, and
C00042.

``` r
library(volcalc)
```

``` r
out_path <- tempdir()
# download a .mol file from KEGG
files <- get_mol_kegg(c("C16181", "C00042"), dir = out_path)
calc_vol(files$mol_path)
#> # A tibble: 2 × 5
#>   mol_path                                          formula name    rvi category
#>   <chr>                                             <chr>   <chr> <dbl> <fct>   
#> 1 /var/folders/wr/by_lst2d2fngf67mknmgf4340000gn/T… C6H7Cl… beta…  6.98 high    
#> 2 /var/folders/wr/by_lst2d2fngf67mknmgf4340000gn/T… C4H6O4  Succ…  2.57 high

#alternatively, supply a SMILES representation
calc_vol(c("C1(C(C(C(C(C1Cl)Cl)Cl)Cl)Cl)O",  "C(CC(=O)O)C(=O)O"), from = "smiles")
#> # A tibble: 2 × 5
#>   smiles                        formula  name    rvi category
#>   <chr>                         <chr>    <chr> <dbl> <fct>   
#> 1 C1(C(C(C(C(C1Cl)Cl)Cl)Cl)Cl)O C6H7Cl5O <NA>   6.98 high    
#> 2 C(CC(=O)O)C(=O)O              C4H6O4   <NA>   2.57 high
```

This returns a dataframe with columns specifying general info about the
compound, and the compound’s calculated volatility and corresponding
volatility category. The functional group counts underlying the
volatility can be additionally returned with `return_fx_groups = TRUE`,
and the intermediate calculation steps with `return_calc_steps = TRUE`.

## Code of Conduct

Please note that the volcalc project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## How to contribute

We appreciate many kinds of feedback and contributions to this R
package. If you find a bug, are interested in an additional feature, or
have made improvements to the package that you want to share, feel free
to file an [issue](https://github.com/Meredith-Lab/volcalc/issues/new)
in this GitHub repo.

## How to cite

If you use this package in your published work, please cite it using the
reference below:

``` r
citation("volcalc")
#> To cite volcalc in publications please use:
#> 
#>   Riemer K, Scott E (2023). _volcalc: Calculate Volatility of Chemical
#>   Compounds_. doi:10.5281/zenodo.8015155
#>   <https://doi.org/10.5281/zenodo.8015155>, R package version 2.0.0.
#> 
#> Please also cite the related manuscript:
#> 
#>   Meredith L, Ledford S, Riemer K, Geffre P, Graves K, Honeker L,
#>   LeBauer D, Tfaily M, Krechmer J (2023). "Automating methods for
#>   estimating metabolite volatility." _Frontiers in Microbiology_.
#>   doi:10.3389/fmicb.2023.1267234
#>   <https://doi.org/10.3389/fmicb.2023.1267234>.
#> 
#> To see these entries in BibTeX format, use 'print(<citation>,
#> bibtex=TRUE)', 'toBibtex(.)', or set
#> 'options(citation.bibtex.max=999)'.
```

### References

Pankow, J.F., Asher, W.E., 2008. SIMPOL.1: a simple group contribution
method for predicting vapor pressures and enthalpies of vaporization of
multifunctional organic compounds. Atmos. Chem. Phys.
<https://doi.org/10.5194/acp-8-2773-2008>
