\name{dbitri}
\alias{dbitri}
\title{Density function for a bi-triangular random variable.}
\usage{
  dbitri(x, a = log2(1.2), b = log2(4), m = log2(2))
}
\arguments{
  \item{x}{vector}

  \item{a}{location of point ... Default a = log2(1.2).}

  \item{b}{location of point ... Default b = log2(4).}

  \item{m}{location of the midpoint of the triangle.
  Default m = log2(2).}
}
\value{
  Gives the density function.
}
\description{
  Density function for a bi-triangular random variable.
}
\details{
  For more details see M. Langaas et al. JRSS B 2005.
}
\examples{
curve(dbitri, -4, 4)
}
\author{
  Maarten van Iterson
}

