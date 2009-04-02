\name{Power}
\alias{Power}

\title{Estimates the average power of the pilot data}

\description{The function \code{Power} estimates the average power of the pilot data at a given false discovery rate. The average power can be estimated for sample sizes other then the pilot data.}

\usage{Power(x, threshold = 0, fdr = 0.1, samplesizes = NULL, 
       plot = FALSE, type = "l", ylim = c(0, 1), xlim = c(0, 1), xlab = "p-value", 
       ylab = "average power", main, sub, \dots)
}

\arguments{
  \item{x}{object of class \code{\link{SampleSize-class}}}
	\item{threshold}{threshold for truncation of the density of effect-sizes. A small symmetrical region around zero can be defined that will be excluded from the density of effect sizes.}
  \item{fdr}{numeric given the desired false discovery rate threshold, possibly a vector like: \code{c(0.1, 0.05)}}
  \item{samplesizes}{numeric vector, giving the samplesizes for which the power will be estimated}
  \item{plot}{logical if \code{TRUE} plots the power curve and intersection for the given fdr(default \code{plot=FALSE}).}
  \item{type}{what type of plot should be drawn}
  \item{ylim}{range of y values}
  \item{xlim}{range of x values}
  \item{xlab}{a title for the x axis}
  \item{ylab}{a title for the y axis}
  \item{main}{an overall title for the plot}
  \item{sub}{a sub title for the plot}
  \item{\dots}{additional arguments to \code{\link{plot}} or \code{\link{par}}}
}

\details{
Having estimated the proportion of non-differentially expressed genes and the density of effect-sizes an estimate of the average power is calculated at a given false discovery rate and sample size. Power and sample size analysis can be performed by estimating the average power for sample sizes other then the pilot data. Ferreira and Zwinderman (2006) show that the average power (the proportion of rejected hypothesis out of all alternative hypothesis) of the adaptive Benjamini-Hochberg procedure is estimated by the point where the CDF of p-values under the alternative hypothesis intersects the straight line with a slope determined by the given FDR and the estimated proportion of non-differentially expressed genes. Using \code{plot=TRUE} the CDF of p-values under the alternative hypothesis is shown with the intersecting lines for different FDRs and sample sizes. Don't use more than 2 FDRs and more than 3 different sample sizes because the plot will then be unclear.
}

\value{
	estimated average power, probabily a named matrix with on the rows different sample sizes and on the columns the different false discovery rates.
}

\references{
Ferreira, F.A., Zwinderman, A., (2006).
Approximate Power and Sample Size Calculations with Microarray Data: An Illustration.
\emph{Statistical Applications in Genetics and Molecular Biology} 5, (1).}

\author{Maarten van Iterson}

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
Power(ss)

##creating a plot estimate power vs. sample size
pwr <- Power(ss, plot = FALSE, samplesizes = c(5, 10, 15, 20), fdr=0.01)
plot(c(5, 10, 15, 20), pwr, ylim = c(0, 1), type = "b", ylab = "Power", xlab = "Sample size per group")
legend("bottomright", colnames(pwr), col=c(1:ncol(pwr)), pch=1, lty=1)

##creating a plot estimate power vs. sample size for different false discovery rates
pwr <- Power(ss, plot = FALSE, samplesizes = c(5, 10, 15, 20), fdr=c(0.01, 0.05))
matplot(c(5, 10, 15, 20), pwr, ylim = c(0, 1), type = "b", pch=1, ylab = "Power", xlab = "Sample size per group")
legend("bottomright", colnames(pwr), col=c(1:ncol(pwr)), pch=1, lty=1)

}

\keyword{hplot}

