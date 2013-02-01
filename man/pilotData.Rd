\name{pilotData}
\alias{pilotData}
\title{User friendly interface to class "PilotData"}
\usage{
  pilotData(statistics = NULL, samplesize = NULL,
    distribution = c("norm", "t", "f", "chisq"), ...)
}
\arguments{
  \item{statistics}{vector of test statistics}

  \item{samplesize}{total sample size of the pilot-data or
  effective sample size in two-group case (see Details for
  more information).}

  \item{distribution}{type of the null/alternative
  distribution, one of 'norm', 't', 'f' or 'chisq'}

  \item{...}{additional arguments for the distribution like
  degrees of freedom}
}
\value{
  object of class "PilotData"
}
\description{
  User friendly interface to class "PilotData"
}
\details{
  details follow In the two-group case the effective sample
  size is defined as the square-root of the inverse of 1/n1
  + 1/n2.
}
\examples{
pd <- pilotData(statistics=rnorm(100), samplesize=10, distribution="norm")
pd
plot(pd)
}
\author{
  Maarten van Iterson
}

