% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/query.R
\name{queryCheaWeb}
\alias{queryCheaWeb}
\title{Web query ChEA3}
\usage{
queryCheaWeb(geneset, set_name = "untitled", n_results = 10,
  libnames = names(chea3::libs), integrate_libs = T, method = c("FET"),
  background = NULL)
}
\arguments{
\item{geneset}{Character array of gene set to query}

\item{set_name}{Name of query. Defaults to "untitled".}

\item{n_results}{Number of results to return for each TF regulon library. Defaults to 10. Input "all"
to return results for all transcription factors.}

\item{libnames}{Character array of TF regulon library names to use. Defaults to all libraries available in the package.}

\item{integrate_libs}{Integrate library results. Often results in a better transcription factor ranking. Defaults to TRUE.}

\item{method}{Method(s) for comparison of TF regulons and query. Defaults to Fisher's Exact Test}

\item{background}{The size of the universe of genes under consideration. User may set to a constant. Default is NULL which results in background being set to the size of the union of the unique genes of the query and TF regulon library.}
}
\value{
A json of datasets with the top TF regulon hits for each TF regulon library tested.
}
\description{
Serves as a wrapper function for queryChea()
}
\examples{
data(example_query)
results = queryChea3(set_name = "exp1_DEG",example_query)

}
