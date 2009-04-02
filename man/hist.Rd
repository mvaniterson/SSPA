\name{hist}
\docType{methods}
\alias{hist}
\alias{hist,PilotData-method}

\title{Histogram of p-values of the pilot data}

\description{The function \code{hist} computes a histogram of p-values of the given pilot data. Two-sided p-values are computed either using the normal or Student t distribution based on the test statistics that are stored in an object of \code{\linkS4class{PilotData}}.}

\usage{
hist(x, \dots)
}

\arguments{
  \item{x}{object of \code{\linkS4class{PilotData}}}
  \item{\dots}{additional arguments to \code{hist}}
}

\details{p-values calculated under the null hypothesis (non-differentially expressed genes) are assumed to be uniformly distributed on [0,1]. p-values calculated under the alternative hypothesis (differentially expressed genes) are assumed to accummelated near zero. The height of the flat part of the histogram is an indication of how many hypothesis are calculated under the null hypothesis.
}

\value{an object of class "histogram", see \code{\link[graphics]{hist}}.}

\author{Maarten van Iterson}

\seealso{\code{\link[graphics]{hist}} and \code{\link{plot}}}

\examples{
library(multtest)
data(golub)
teststat <- mt.teststat(golub, golub.cl)
table(golub.cl)
pd <- pilotData(name="golub", testStatistics=teststat, sampleSizeA=11, sampleSizeB=27)
hist(pd)
}

\keyword{hplot}

