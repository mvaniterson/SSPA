\name{PilotData-class}
\docType{class}
\alias{PilotData-class}
\alias{distribution,PilotData-method}
\alias{dof,PilotData-method}
\alias{effectivesamplesize,PilotData-method}
\alias{plot,PilotData-method}
\alias{pvalues,PilotData-method}
\alias{show,PilotData-method}
\alias{statistics,PilotData-method}

\title{Class PilotData}
\description{The \code{"\linkS4class{PilotData}"} class is a container class. It contains all necessary information of the pilot data for performing the power and sample size analysis.}
\section{Objects from the Class}{
Objects can be created by calling the function \code{\link{pilotData}}.
}
\section{Slots}{
	 \describe{    
    \item{\code{statistics}:}{vector of test statistics, type \code{"numeric"}}
    \item{\code{effectivesamplesize}:}{in the two group case effective sample size \eqn{(\frac{1}{n_A}+\frac{1}{n_B})^{-1}} with \eqn{n_A} the sample size in group A, type \code{"numeric"}}
    \item{\code{pvalues}:}{vector of p-values, type \code{"numeric"}}
    \item{\code{dof}:}{degrees of freedom only used for Student t distribution, type \code{"numeric"}}
    \item{\code{distribution}:}{null distribution either Normal or Student t, type \code{"character"}}
  }
}
\section{Methods}{
  \describe{
    \item{distribution}{\code{signature(object = "PilotData")}: ... }
    \item{dof}{\code{signature(object = "PilotData")}: ... }
    \item{effectivesamplesize}{\code{signature(object = "PilotData")}: ... }
    \item{plot}{\code{signature(x = "PilotData")}: ... }
    \item{pvalues}{\code{signature(object = "PilotData")}: ... }
    \item{show}{\code{signature(object = "PilotData")}: ... }
    \item{statistics}{\code{signature(object = "PilotData")}: ... }
	 }
}

\author{Maarten van Iterson}

\seealso{
\code{\link{pilotData}},
\code{\link{SampleSize-class}}}

\keyword{classes}
