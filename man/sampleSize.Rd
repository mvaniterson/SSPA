\name{sampleSize}
\alias{sampleSize}

\title{Creates an object of class SampleSize}

\description{The function \code{sampleSize} initializes a \code{\linkS4class{SampleSize}}-object. The density of effect-sizes is estimated using a deconvolution estimator or constrained optimization using B-splines. 
}

\usage{
sampleSize(PilotData, method = c("Langaas", "Storey", "Ferreira", "Ruppert", "Userdefined"), from = -6, to = 6, resolution = 2^10, kernel = c("fan", "wand", "sinc"), pi0 = seq(0.1, 0.99, 0.01), adjust=TRUE, a=0.5, bandwidth, nKnots = 11, bDegree = 3, ...)
}

\arguments{
  \item{PilotData}{object of class \code{\linkS4class{PilotData}} on which the sample size and power analysis will be performed}
  \item{method}{character string giving the method for estimation of the fraction of non-differentially expressed genes either one of:
		\itemize{
			\item "Langaas" (default) 
      \item "Storey" 
      \item "Ferreira" 
      \item "Ruppert"  
      \item "Userdefined"}}
  \item{from}{Lower bound of the density of effect-sizes (the range should be symmetric)}
  \item{to}{Upper bound of the density of effect-sizes.}
  \item{resolution}{the number of points on which the density of effect-sizes will be estimated (must be a power of 2)}
  \item{kernel}{the kernel type used in the deconvolution estimator either one of:
		\itemize{
    	\item "fan" (default)
      \item "wand"
      \item "sinc"}}
  \item{pi0}{numeric or a vector of type numeric giving the fraction of non-differentially expressed genes. If \code{method="Userdefined"} a userdefined pi0 is obligated. If method="Ferreira" a range of values should be given e.g. \code{seq(0.01, 0.99, 0.01)}}
	\item{adjust}{is the density of effect sizes adjusted (\code{adjust=TRUE}, default) or not  (\code{adjust=FALSE})}
	\item{a}{For better upper bound estimates of pi0 as described by Efron \emph{et al.} (2001).}
	\item{bandwidth}{Bandwidth for the kernel estimator, default is the plugging bandwidth of Fan.}
	\item{nKnots}{ number of knots used in Rupperts method for estimating the proportion of non-differentially expressed genes and density of effect-sizes}
	\item{bDegree}{degree of B-spline basis used in Rupperts method for estimating the proportion of non-differentially expressed genes and density of effect-sizes}
  \item{\dots}{Additional parameters for method for the estimation of the fraction of non-differentially expressed genes (Doesn't work yet!)}
}

\details{ 
The \code{sampleSize} function performes the estimation of the proportion of non-differentially expressed genes using one of the three methods, \code{"Langaas", "Storey", "Ferreira"} if \code{method="Userdefined"} a userdefined \code{pi0} is needed and estimation of the proportion of non-differentially expressed genes will be skipped. A deconvolution estimator is implemented using the Fast Fourier Transform Algorithm \code{fft()} for estimations of the density of effect sizes. 
If \code{method="Ruppert"} constrained optimization using B-splines is used, for this method two additional packages needs to be installed namely \code{quadprog} and \code{splines}. Both the proportion of non-differentially expressed genes and the density of effect sizes are estimated with Ruppert's method. In contrast to the orginal method of Ruppert et al. (2007) we made a modification on the estimation of the density of effect sizes allowing for negative effect sizes as well.
}

\value{Object of class \code{\linkS4class{SampleSize}}}

\references{
Ferreira, F.A., Zwinderman, A., (2006).
Approximate Power and Sample Size Calculations with Microarray Data: An Illustration.
\emph{Statistical Applications in Genetics and Molecular Biology} 5, (1).

Ferkingstad, E., Langaas, M., and Lindqvist, B. (2005).
Estimating the proportion of true null hypotheses, with application
to DNA microarray data.
\emph{Journal of the Royal Statistical Society Series} B, 67, 555-572.

Storey, J.D., (2002).
A direct approach to false discovery rates.
\emph{J.R. Statist. Med.} 27, 1960-1972.

Ruppert, D. and Nettleton, D. and Hwang, J.T.G., (2007).
Exploring the information in p-values for the analysis and planning of multiple-test experiments.
\emph{Biometrics}, 63, 2, 483-95.

Efron, B.  and Tibshirani, R. and Storey, J.D. and Tusher, V. (2001)
Empirical Bayes Analysis of a Microarray Experiment.
\emph{Journal of the American Statistical Association}, 96, 456, 1151-1160.
}

\author{Maarten van Iterson}

\seealso{\code{\link[limma]{convest}}, \code{\link[qvalue]{qvalue}} and \code{\linkS4class{SampleSize}} use \code{class?SampleSize}}

\examples{
library(multtest)
data(golub)
teststat <- mt.teststat(golub, golub.cl)
table(golub.cl)
X <- model.matrix(~golub.cl)
effectivesamplesize <- 1/solve(crossprod(X))[2,2]
pd <- pilotData(statistics=teststat, effectivesamplesize=effectivesamplesize)
plot(pd)
ss <- sampleSize(pd)
}

\keyword{classes}

