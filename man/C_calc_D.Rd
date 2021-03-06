% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{C_calc_D}
\alias{C_calc_D}
\title{Computes mean density for each occasion}
\usage{
C_calc_D(D, J, pr0, tpms)
}
\arguments{
\item{D}{density parameter}

\item{J}{total number of occasions}

\item{pr0}{initial distribution over life states}

\item{tpms}{list of transition probability matrices (one per occasion)}
}
\value{
vector with j^th entry the mean density on occasion j
}
\description{
Computes mean density for each occasion
}
