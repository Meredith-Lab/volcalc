## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

The note on check is:

>Package has a FOSS license but eventually depends on the following
>  package which may restrict use:
>    ChemmineOB

ChemmineOB is an R package with the [Artistic-2.0 license](https://github.com/girke-lab/ChemmineOB/blob/master/LICENSE), which does appear to be FOSS.  Other R packages

There was an error on the CRAN pre-test on an initial submission on Debian:

>* checking package dependencies ... ERROR
Package required but not available: ‘ChemmineOB’

I've confirmed that `ChemmineOB` is available and can be installed on Debian using rocker/r-devel. I assume this actually means that the installation of `ChemmineOB` errored due to missing system dependencies (OpenBabel and eigen3).  In response, we've moved `ChemmineOB` from Imports to Suggests, checked that it is installed before functions from it are called, and conditionally skipped tests and examples when `ChemmineOB` is not installed.
