## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

The note on check is:

>Package has a FOSS license but eventually depends on the following
>  package which may restrict use:
>    ChemmineOB

ChemmineOB is an R package with the [Artistic-2.0 license](https://github.com/girke-lab/ChemmineOB/blob/master/LICENSE), which does appear to be FOSS.  Other R packages

Another note came up in a CRAN submission pre-test only on Debian

>* checking package dependencies ... ERROR
Package required but not available: ‘ChemmineOB’

I've confirmed that `ChemmineOB` can be installed on Debian using rocker/r-devel. It is possible that Bioconductor might have been temporarily unavailable at the moment this check was run.  It is also possible that there was an error installing ChemmineOB because its SystemRequirements (OpenBabel and eigen3) were not satisfied on the runner.  I'm hoping its former, but if this pre-test fails again, I'll look into fixes (likely moving ChemmineOB into Suggests).
