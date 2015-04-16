\name{simdat}
\alias{simdat}
\title{Generate simulated microarray data using the bitriangular distribution.}
\usage{
  simdat(mu, m, pi0, J, nullX = function(x) rnorm(x, 0, 1),
    nullY = function(x) rnorm(x, 0, 1), noise = 0.01)
}
\arguments{
  \item{mu}{vector of effect sizes drawn from the
  bitriangular distribution.}

  \item{m}{number of features (genes, tags, ...).}

  \item{pi0}{proportion of nondifferentially expressed
  features.}

  \item{J}{number of samples per group.}

  \item{nullX}{the distribution of nondifferentially
  expressed features.}

  \item{nullY}{the distribution of nondifferentially
  expressed features.}

  \item{noise}{standard deviation of the additive noise.}
}
\value{
  Matrix of size m x (2J), containing the simulated values.
}
\description{
  Simulated microarray data.
}
\details{
  details follow
}
\examples{
##generate two-group microarray data
m <- 5000 ##number of genes
J <- 10 ##sample size per group
pi0 <- 0.8 ##proportion of non-differentially expressed genes
m0 <- as.integer(m*pi0)
mu <- rbitri(m - m0, a = log2(1.2), b = log2(4), m = log2(2)) #effect size distribution
data <- simdat(mu, m=m, pi0=pi0, J=J, noise=0.01)
}
\author{
  Maarten van Iterson
}

