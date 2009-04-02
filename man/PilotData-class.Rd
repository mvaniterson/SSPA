\name{PilotData-class}
\docType{class}
\alias{PilotData-class}

\title{Class PilotData}
\description{The \code{"\linkS4class{PilotData}"} class is a container class. It contains all necessary information of the pilot data for performing the power and sample size analysis.}
\section{Objects from the Class}{
Objects can be created by calling the function \code{\link{pilotData}}.
}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{name of the pilot data, type \code{"character"}}
    \item{\code{testStatistics}:}{vector of test statistics, type \code{"numeric"}}
    \item{\code{sampleSize}:}{effective sample size \eqn{(\frac{1}{n_A}+\frac{1}{n_B})^{-1}} with \eqn{n_A} the sample size in group A, type \code{"numeric"}}
    \item{\code{pValues}:}{vector of p-values, type \code{"numeric"}}
    \item{\code{dof}:}{degrees of freedom only used for Student t distribution, type \code{"numeric"}}
    \item{\code{nullDist}:}{null distribution either normal or Student t, type \code{"character"}}
  }
}
\section{Methods}{
\code{\link{show,PilotData-method}}, \code{\link{hist}} and \code{\link{plot}}
}

\author{Maarten van Iterson}

\seealso{
\code{\link{pilotData}},
\code{\link{SampleSize-class}}}

\keyword{classes}
