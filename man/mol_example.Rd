% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mol_example.R
\name{mol_example}
\alias{mol_example}
\title{Example .mol files}
\usage{
mol_example()
}
\value{
File paths to installed example .mol files.
}
\description{
\code{volcalc} comes bundled with some example .mol files in its \code{inst/extdata}
directory.  This function provides easy access to them.
}
\details{
File names are the KEGG identifiers.  Compound names are as follows:
\itemize{
\item C00031: D-Glucose
\item C00157: Phosphatidylcholine
\item C08491: (-)-Jasmonic acid
\item C16181: beta-2,3,4,5,6-Pentachlorocyclohexanol
\item C16286: Geosmin
\item C16521: Isoprene
}
}
\examples{
#return paths to all example .mol files
mol_example()

#examine the contents of a file
readLines(mol_example()[1])
}
