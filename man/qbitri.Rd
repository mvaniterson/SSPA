\name{qbitri}
\alias{qbitri}
\title{Quantile function for a bi-triangular random variable.}
\usage{
  qbitri(p, a = log2(1.2), b = log2(4), m = log2(2))
}
\arguments{
  \item{p}{vector of probabilities.}

  \item{a}{location of point, ... Default a = log2(1.2).}

  \item{b}{location of point, ... Default b = log2(4).}

  \item{m}{location of the midpoint of the triangle.
  Default m = log2(2).}
}
\value{
  Gives the quantile function.
}
\description{
  Quantile function for a bi-triangular random variable.
}
\details{
  For more details see M. Langaas et al. JRSS B 2005.
}
\examples{
curve(qbitri, 0, 1)
}
\author{
  Maarten van Iterson
}

