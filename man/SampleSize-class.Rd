\name{SampleSize-class}
\docType{class}
\alias{SampleSize-class}
\alias{lambda,SampleSize-method}
\alias{pi0,SampleSize-method}
\alias{plot,SampleSize-method}
\alias{show,SampleSize-method}
\alias{theta,SampleSize-method}

\title{Class SampleSize}
\description{By contructing an object of \code{SampleSize} the power and sample size analysis will be performed.}
\section{Objects from the Class}{
Objects can be created by calling the function \code{\link{sampleSize}}.

}
\section{Slots}{
	 \describe{
    \item{\code{pi0}:}{proportion of non-differentially expressed genes, \code{"named-list"}}
    \item{\code{lambda}:}{vector containing the estimated density of effect sizes, \code{"numeric"}}
    \item{\code{theta}:}{vector containing the range over which the density of effect sizes is estimated, \code{"numeric"}}
    \item{\code{adjust}:}{is the density of effect sizes adjusted (\code{adjust=TRUE}) or not  (\code{adjust=FALSE}), \code{"logical"}}
    \item{\code{method}:}{method used for the estimation of the density of effect sizes, \code{"character"}}
    \item{\code{bandwidth}:}{Default bandwidth given by \eqn{log_{10}(n)^{-\frac{1}{2}}} where \eqn{n} is the number of genes.}
    \item{\code{kernel}:}{kernel used for the deconvolution estimator, \code{"character"}}
    \item{\code{nKnots}:}{number of knots used for the B-splines (only visible if \code{method="Ruppert"})}
    \item{\code{bDegree}:}{B-spline degree(only visible if \code{method="Ruppert"})}
    \item{\code{statistics}:}{vector of test statistics, type \code{"numeric"}}
    \item{\code{effectivesampleSize}:}{effective sample size \eqn{(\frac{1}{n_A}+\frac{1}{n_B})^{-1}} with \eqn{n_A} the sample size in group A, type \code{"numeric"}}
    \item{\code{pvalues}:}{vector of p-values, type \code{"numeric"}}
    \item{\code{dof}:}{degrees of freedom only used for Student t distribution, type \code{"numeric"}}
    \item{\code{distribution}:}{null distribution either normal or Student t, type \code{"character"}}
  }
}
\section{Extends}{
Class \code{\linkS4class{PilotData}}, directly.
}
\section{Methods}{
  \describe{
    \item{lambda}{\code{signature(object = "SampleSize")}: ... }
    \item{pi0}{\code{signature(object = "SampleSize")}: ... }
    \item{plot}{\code{signature(x = "SampleSize")}: ... }
    \item{show}{\code{signature(object = "SampleSize")}: ... }
    \item{theta}{\code{signature(object = "SampleSize")}: ... }
	 }
}

\author{Maarten van Iterson}

\seealso{
	\code{\link{pilotData}},
	\code{\linkS4class{PilotData}}
}
\keyword{classes}
