# volcalc (development version)

* Change output of `get_fx_groups()` and `calc_vol()` from data frame to tibble
* `get_fx_groups()` and `calc_vol()` no longer depend on KEGG or take KEGG compound IDs or pathway IDs.  Instead, `calc_vol()` accepts a path to a .mol file as input.
* `calc_vol()` is vectorized and accepts multiple compounds as input.
* Moved SIMPOL.1 calculations out of `calc_vol()` and into to their own function, `simpol1()`, to pave the way for future expansions using other methods.  The "manual" workflow is now .mol file |> `ChemmineR::read.SDFset()` |> `get_fx_groups()` |> `simpol1()`
* The output of `calc_vol()` (and `simpol1()`) now contains a column called `log10_P` instead of `log_Sum`, equivalent to `log_Sum` + the coefficient for b_0(T)

# volcalc 1.0.2

* Minor change in calculation in `calc_vol()`---remove amines functional group to avoid double counting with primary amines (#49)

# volcalc 1.0.1

* Minor change in calculation in `calc_vol()`---use 293.15K for temperature (instead of 293K) to match the temperature used in Pankow & Asher (2008)

# volcalc 1.0.0

* Initial release of `volcalc`.  This is the version of the code used in an in-prep manuscript.  Version 2.0.0 will include **breaking changes**

# volcalc 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
