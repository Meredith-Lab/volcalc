# volcalc (development version)

* Added pkgdown website
* `get_fx_groups()` now only counts the smallest set of smallest rings (#57)
* Fixed a bug that caused the number of aromatic hydroxyl groups (e.g. phenols) to be miscounted (#46)

# volcalc 1.0.2


* Minor change in calculation in `calc_vol()`---remove amines functional group to avoid double counting with primary amines (#49)

# volcalc 1.0.1

* Minor change in calculation in `calc_vol()`---use 293.15K for temperature (instead of 293K) to match the temperature used in Pankow & Asher (2008)

# volcalc 1.0.0

* Initial release of `volcalc`.  This is the version of the code used in an in-prep manuscript.  Version 2.0.0 will include **breaking changes**

# volcalc 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
