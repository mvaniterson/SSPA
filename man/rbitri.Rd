\name{rbitri}
\alias{rbitri}
\title{Random generation of bitriangular distributed values.}
\usage{
  rbitri(n, a = log2(1.2), b = log2(4), m = log2(2))
}
\arguments{
  \item{n}{number of observations.}

  \item{a}{location of point, ... Default a = log2(1.2).}

  \item{b}{location of point, ... Default b = log2(4).}

  \item{m}{location of the midpoint of the triangle.
  Default m = log2(2).}
}
\value{
  Generates random deviates.
}
\description{
  Random generation of bitriangular distributed values.
}
\details{
  For more details see M. Langaas et al. JRSS B 2005.
}
\examples{
hist(rbitri(100), freq=FALSE)
curve(dbitri, add=TRUE)
}
\author{
  Maarten van Iterson
}

