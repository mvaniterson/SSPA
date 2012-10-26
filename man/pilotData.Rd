\name{pilotData}
\alias{pilotData}

\title{Creates an object of class "PilotData"}

\description{The function \code{pilotData} initializes a \code{\linkS4class{PilotData}}-object. Information of the pilot data and the null distribution is stored e.g. name of pilot experiment, test statistics and number of samples used.}

\usage{pilotData(statistics = NULL, effectivesamplesize = NULL, dof=double(1), distribution = c("Normal", "Student"))}

\arguments{
  \item{statistics}{vector of type numeric containing the set of test statistics}
  \item{effectivesamplesize}{the effective sample size}
	\item{dof}{degrees of freedom, only necessary when \code{distribution=="Student"}}
  \item{distribution}{distribution under the null hypothesis either one of:
    \itemize{
      \item "Normal",
      \item "Student",       
    }}}

\details{Based on the given null distribution two-sided p-values are calculated from the test statistics. Some additional checks on the data are performed. Once an object of \code{\linkS4class{PilotData}} is created the \code{\link{sampleSize}}-function can be called.
}

\value{object of class \code{\linkS4class{PilotData}}}

\author{Maarten van Iterson}

\seealso{\code{plot} and \code{\linkS4class{PilotData}} use \code{class?PilotData}}

\examples{
library(multtest)
data(golub)
teststat <- mt.teststat(golub, golub.cl)
table(golub.cl)
X <- model.matrix(~golub.cl)
effectivesamplesize <- 1/solve(crossprod(X))[2,2]
pd <- pilotData(statistics=teststat, effectivesamplesize=effectivesamplesize)
}

\keyword{classes}
