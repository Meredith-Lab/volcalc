# volcalc (development version)

* It is now possible to supply input to `calc_vol()` as a vector of SMILES strings with `from = "smiles"`
* Users can now choose from RVI thresholds for non-volatile, low, moderate, and high volatility for clean atmosphere, polluted atmosphere, or soil using the `environment` parameter of `calc_vol()
* A coefficient for amides has been removed from the "Meredith" method of `simpol1()` to avoid double-counting amides.

## Changes to `get_mol_kegg()`

* The `pathway_ids` argument of `get_mol_kegg()` now also accepts pathway *module* IDs (e.g. "M00082")
* `get_mol_kegg()` got a significant speed improvement (#84)
* `get_mol_kegg()` will skip downloading a .mol file if it is already present by default (override with `force=TRUE`)

## Changes to `get_fx_groups()`

* `mass` column renamed to `molecular_weight` and addition of an `exact_mass` column
* change to how non-aromatic carbon double bonds are counted.  Now using SMARTS pattern "C=C"
* now returns `hydroxyl_total` and `hydroxyl_aliphatic` instead of `hydroxyl_groups`
* now returns `rings_total` and `rings_aliphatic` instead of `rings`
* counts additional groups: 
  - aromatic amines
  - primary, secondary, and tertiary amides
  - hydroperoxides
  - carbonylperoxyacids
  - carbonylperoxynitrates
  - alkyl, alicyclic, and aromatic ethers (in addition to total ethers)
  - nitrophenols
  - nitroesters
* changed `ether` to `ether_alkyl` and added `ether_total` (matching any R-O-R)

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
