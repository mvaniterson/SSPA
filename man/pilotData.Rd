\name{pilotData}
\alias{pilotData}

\title{Creates an object of class "PilotData"}

\description{The function \code{pilotData} initializes a \code{\linkS4class{PilotData}}-object. Information of the pilot data and the null distribution is stored e.g. name of pilot experiment, test statistics and number of samples used.}

\usage{pilotData(name = "Unknown Experiment", testStatistics = double(1), sampleSizeA = double(1), sampleSizeB = double(1), dof = double(1),	 nullDist = c("normal", "student"))}

\arguments{
  \item{name}{character string giving the experiment name}
  \item{testStatistics}{vector of type numeric containing the set of test statistics}
  \item{sampleSizeA}{the samplesize of group A}
  \item{sampleSizeB}{the samplesize of group B}
  \item{dof}{degree of freedom for a Student t distribution}
  \item{nullDist}{distribution under the null hypothesis either one of:
    \itemize{
      \item "normal",
      \item "student", (only in case of method="Ruppert" see \code{\link{sampleSize}})      
    }}}


\details{Based on the given null distribution two-sided p-values are calculated from the test statistics. Some additional checks on the data are performed. Once an object of \code{\linkS4class{PilotData}} is created the \code{\link{sampleSize}}-function can be called.
}

\value{object of class \code{\linkS4class{PilotData}}}

\author{Maarten van Iterson}

\seealso{\code{hist}, \code{plot} and \code{\linkS4class{PilotData}} use \code{class?PilotData}}

\examples{
library(multtest)
data(golub)
teststat <- mt.teststat(golub, golub.cl)
table(golub.cl)
pd <- pilotData(name="golub", testStatistics=teststat, sampleSizeA=11, sampleSizeB=27)
}

\keyword{classes}
