# volcalc (development version)

* Change output of `get_fx_groups()` and `calc_vol()` from data frame to tibble
* `get_fx_groups()` and `calc_vol()` no longer depend on KEGG or take KEGG compound IDs or pathway IDs.  Instead, `calc_vol()` accepts a path to a .mol file as input.
* `calc_vol()` is vectorized and accepts multiple compounds as input.

# volcalc 1.0.0

* Initial release of `volcalc`.  This is the version of the code used in an in-prep manuscript.  Version 2.0.0 will include **breaking changes**

# volcalc 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
