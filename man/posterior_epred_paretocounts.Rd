% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paretocounts.R
\name{posterior_epred_paretocounts}
\alias{posterior_epred_paretocounts}
\title{Arrange data for posterior epred prediction}
\usage{
posterior_epred_paretocounts(prep)
}
\arguments{
\item{prep}{see documentation from `get_dpar` within `brms`}
}
\value{
NA
}
\description{
This function is not called directly. It is only used to allow add_epred_draws() from the
`tidybayes` package to work with a `brmfit` object with a paretocounts_lpdf.
}
\examples{
NA
}
