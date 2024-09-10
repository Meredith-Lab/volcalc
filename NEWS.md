# volcalc (development version)

* adds a `validate = TRUE` option to `calc_vol()` and `get_fx_groups()` that returns `NA`s when there are suspected errors in parsing SMILES or .mol files. This is unfortunately not available on Windows due to differences in the windows version of `ChemmineOB`
* `KEGGREST` is no longer a dependency of `volcalc` (previously used in `get_mol_kegg()`)

# volcalc 2.1.2

* There is no release for this version as it was a rejected CRAN submission.

# volcalc 2.1.1

* Added Dr. Laura Meredith as a package author as the concept for volcalc originated from her.
* Added package contributors S. Marshall Ledford (our main user and beta tester since day 1!) and Tamás Stirling (for help with SMARTS patterns)
* Updated citation files accordingly

# volcalc 2.1.0

* The manuscript associated with `volcalc` is now published in Frontiers in Microbiology 🎉. DOI: 10.3389/fmicb.2023.1267234
* There are now two vignettes available which can be viewed with `browseVignettes("volcalc")` or on the package website

## Miscelanous changes

* New example .mol files were added.  See `?mol_example()`
* `mol_example()` no longer takes any arguments and just returns file paths to all example .mol files

## Changes to `calc_vol()` and `simpol1()`

* It is now possible to supply input to `calc_vol()` as a vector of SMILES strings with `from = "smiles"`.
* Users can now choose from RVI thresholds for non-volatile, low, moderate, and high volatility for clean atmosphere, polluted atmosphere, or soil using the `environment` parameter of `calc_vol()`.
* A coefficient for amides has been removed from the "Meredith" method of `simpol1()` to avoid double-counting amides.
* The default for the `method` argument to `calc_vol()` has now been renamed to `"meredith"`.  `"simpol1"` now uses the original SIMPOL.1 method without additional coefficients added in Meredith et al. (2023).
* `simpol1()` gains an argument `meredith` that controls whether just the functional groups in the original SIMPOL.1 method (Pankow & Asher, 2008) is used or if additional coefficients used in Meredith et al. (2023) are also included.
* `simpol1()` now takes into account amide functional groups.

## Changes to `get_mol_kegg()`

* The `pathway_ids` argument of `get_mol_kegg()` now also accepts pathway *module* IDs (e.g. "M00082").
* `get_mol_kegg()` got a significant speed improvement (#84).
* `get_mol_kegg()` will skip downloading a .mol file if it is already present by default (override with `force=TRUE`).

## Changes to `get_fx_groups()`

* `mass` column renamed to `molecular_weight` and addition of an `exact_mass` column.
* change to how non-aromatic carbon double bonds are counted.  Now using SMARTS pattern "C=C"
* now returns `hydroxyl_total` and `hydroxyl_aliphatic` instead of `hydroxyl_groups`
* now returns `rings_total` and `rings_aliphatic` instead of `rings`
* Now captures all functional groups in the SIMPOL.1 method except "number of carbons on the acid side of an amide". This version adds these previously missing groups used by `simpol1()`: 
  - C=C-C=O in a non-aromatic ring
  - non-aromatic carbon double bonds
  - aromatic amines
  - primary, secondary, and tertiary amides
  - hydroperoxides
  - carbonylperoxyacids
  - carbonylperoxynitrates
  - alkyl, alicyclic, and aromatic ethers (in addition to total ethers)
  - nitrophenols
  - nitroesters
* changed `ether` to `ether_alkyl` and added `ether_total` (matching any R-O-R)
* slight change to the SMARTS pattern to capture sulfonate groups to also capture conjugate sulfonic acids

# volcalc 2.0.0

This version includes big (breaking) changes in how the package works!  Please
read the changelog below carefully and check function documentation and examples
to see the new usage of functions.

* Change output of `get_fx_groups()` and `calc_vol()` from data frame to tibble
* `get_fx_groups()` and `calc_vol()` no longer depend on KEGG or take KEGG compound IDs or pathway IDs.  Instead, `calc_vol()` accepts a path to a .mol file as input.
* `calc_vol()` is vectorized and accepts multiple compounds as input.
* Moved SIMPOL.1 calculations out of `calc_vol()` and into to their own function, `simpol1()`, to pave the way for future expansions using other methods.  The "manual" workflow is now .mol file |> `ChemmineR::read.SDFset()` |> `get_fx_groups()` |> `simpol1()`
* The output of `calc_vol()` (and `simpol1()`) now contains a column called `log10_P` instead of `log_Sum`, equivalent to `log_Sum` + the coefficient for b_0(T)
* Output of `calc_vol()` now contains a column with the inputs, named whatever is supplied to `from` (eg. a column called `mol_path` containing paths to mol files)
* A new function, `get_mol_kegg()`, replaces `save_compound_mol()` for downloading mol files from KEGG
* Added pkgdown website
* `get_fx_groups()` now only counts the smallest set of smallest rings (#57)
* Fixed a bug that caused the number of phenols to be miscounted. Rather than counting phenols, `get_fx_groups` now counts aromatic hydroxyl groups (e.g. phenols) to more closely align with Pankow & Asher (2008) (#46)
* package now has a hex logo!
* Fixes a bug in `volcalc` introduced by a bug-fix in `ChemmineR` v3.53.1 (#54)
* The `volatility` column in the output of `calc_vol()` has been renamed to `rvi` (relative volatility index)

# volcalc 1.0.2

* Minor change in calculation in `calc_vol()`---remove amines functional group to avoid double counting with primary amines (#49)

# volcalc 1.0.1

* Minor change in calculation in `calc_vol()`---use 293.15K for temperature (instead of 293K) to match the temperature used in Pankow & Asher (2008)

# volcalc 1.0.0

* Initial release of `volcalc`.  This is the version of the code used in an in-prep manuscript.  Version 2.0.0 will include **breaking changes**

# volcalc 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
