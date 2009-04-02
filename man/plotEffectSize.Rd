\name{plotEffectSize}
\alias{plotEffectSize}

\title{Plots the density of effect sizes of the pilot data}

\description{The function \code{plotEffectSize} plots density of effect sizes of the pilot data.}

\usage{plotEffectSize(x, threshold = 0, xlab = "effect size", ylab = "density of effect sizes", main, sub, \dots)}

\arguments{
  \item{x}{object of class \code{\link{SampleSize-class}}}
	\item{threshold}{threshold for truncation of the density of effect-sizes. The threshold will be taken symmetrical around the y-axis.}
  \item{xlab}{a title for the x axis}
  \item{ylab}{a title for the y axis}
  \item{main}{an overall title for the plot}
  \item{sub}{a sub title for the plot}
  \item{\dots}{additional arguments given to \code{\link{plot}} or \code{\link{par}}}
}

\details{The density of effect sizes describes the effects observed in the pilot data. Usually a bimodal density is observed representing up- and down-regulated genes. The way in which the test statistics is calculated determines what is meant by up- and down-regulation. A small symmetrical region around zero can be defined that will be excluded from the density of effect sizes and thereby increases the estimated average power. }

\references{
Ferreira, F.A., Zwinderman, A., (2006).
Approximate Power and Sample Size Calculations with Microarray Data: An Illustration.
\emph{Statistical Applications in Genetics and Molecular Biology} 5, (1).}

\examples{
library(multtest)
data(golub)
teststat <- mt.teststat(golub, golub.cl)
table(golub.cl)
pd <- pilotData(name="golub", testStatistics=teststat, sampleSizeA=11, sampleSizeB=27)
hist(pd)
plot(pd)
ss <- sampleSize(pd)
plotEffectSize(ss)
}

\keyword{hplot}

