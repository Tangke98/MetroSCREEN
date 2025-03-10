% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kpc.r
\name{regrVonPS}
\alias{regrVonPS}
\title{Check if variable can be regressed to independence on its parents}
\usage{
regrVonPS(G, V, S, suffStat, indepTest = kernelCItest, alpha = 0.2)
}
\arguments{
\item{G}{adjacency matrix, for the graph}

\item{V}{integer, node which we regress}

\item{S}{integer(s), set we regress on}

\item{suffStat}{sufficient statistics to perform the independence test \link{kernelCItest}}

\item{indepTest}{independence test to check for dependence between residuals of V and S}

\item{alpha}{numeric cutoff for significance level of individual partial correlation tests}
}
\value{
regrVonPS() returns the number of p-values smaller than the cutoff, i.e 0 means residuals of V are independent of all variables in S
}
\description{
Uses the generalised additive model \link{gam} to non-linearly and non-parametrically regress variable V on its parents and set of variables S.
}
\author{
Petras Verbyla (\email{petras.verbyla@mrc-bsu.cam.ac.uk})
}
