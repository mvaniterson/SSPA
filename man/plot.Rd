\name{plot}
\docType{methods}
\alias{plot}
\alias{plot-methods}
\alias{plot,ANY,ANY-method}
\alias{plot,PilotData,missing-method}

\title{Empirical cumulative distribution function of p-values of the pilot data}

\description{The function \code{plot} plots an empirical cumulative distribution function of p-values of the pilot data. The test statistics that are given as input to \code{\link{PilotData-class}} are transformed to two-sided p-values either using the normal or Student t distribution.}

\usage{plot(x, y, \dots)}

\arguments{
  \item{x}{object of class \code{\link{PilotData-class}}}
  \item{y}{unused item}
  \item{\dots}{additional arguments given to \code{\link[graphics]{plot}} or \code{\link{par}}}
}

\details{
Empirical cumulative distribution function of p-values. The line at angle of 45 degrees represents the theoretical CDF of a uniform distribution as
expected when all genes are non-differentially expressed. An accumulation of p-values near zero indicates a certain number of differentially expressed genes.
}

\author{Maarten van Iterson}

\examples{
library(multtest)
data(golub)
teststat <- mt.teststat(golub, golub.cl)
table(golub.cl)
pd <- pilotData(name="golub", testStatistics=teststat, sampleSizeA=11, sampleSizeB=27)
hist(pd)
plot(pd)
}

\keyword{hplot}

