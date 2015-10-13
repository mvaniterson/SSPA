\name{sampleSize}
\alias{sampleSize}
\title{User friendly interface to class 'SampleSize'}
\usage{
  sampleSize(PilotData,
    method = c("deconv", "congrad", "tikhonov", "ferreira"),
    control = list(from = -6, to = 6, resolution = 2^9))
}
\arguments{
  \item{PilotData}{object of class 'PilotData'.}

  \item{method}{estimation method one of 'deconv',
    'congrad', 'tikhonov' or 'ferreira'. See 'Details'.}

  \item{control}{A list of control parameters. See
    'Details'.}
}
\value{
  object of class SampleSize.
}
\description{
  User friendly interface to class "SampleSize"
}
\details{
  The default method is 'deconv' which is an kernel deconvolution density estimator implementated using \code{\link{fft}}.
  The 'nncg' is a nonnegative conjugate gradient algorithm based on R's
  implementation see \code{\link{optim}}.
  'tikonov' implements ridge-regression with optimal penalty selection using the L-curve approach. Higher order
  penalties are possible as well using a transformation to standard form
  (see Hansen).

  The 'control' argument is a list that can supply any of the following components. Per method logical checks are performed.

  \itemize{
    \item{deconv:}{
      \itemize{
        \item{method:}{'deconv', 'ferreira'}
        \item{pi0Method:}{the pi0 estimation method one of 'Langaas', 'Storey', 'Ferreira', 'Userdefined'}
        \item{pi0:}{if method = 'ferreira' grid pi0-value need to be suppled e.g. seq(0.1, 0.99, 0.01)}
        \item{adjust:}{Default TRUE, adjust pi0 esitmate if density of effect size is somewhere negative.}
        \item{a:}{Adjust pi0 better approach suggested by Efron. Symmetric range around zero of size 0.5.}
        \item{bandwith:}{Default NULL uses 1/sqrt(log(length(statistics)))}
        \item{kernel:}{Either 'fan', 'wand', 'sinc' kernels can be used.}
        \item{from:}{Density of effect sizes should be estimated from = -6}
        \item{to:}{ to = 6}
        \item{resolution:}{Density of effect sizes should be estimated on 2^9 points.}
        \item{verbose:}{Default FALSE if TRUE additional information is printed to the console.}}
    }
    \item{congrad:}{
      \itemize{
        \item{integration:}{'midpoint', 'trapezoidal', 'simpson'}
        \item{scale:}{'pdfstat', 'cdfstat', 'cdfpval'}
        \item{trim:}{0.01, 0.99}
        \item{symmetric:}{TRUE}
        \item{bin:}{'epdf', 'ecdf'}
        \item{from:}{-6}
        \item{to:}{6}
        \item{resolution:}{500}
        \item{verbose:}{Default FALSE if TRUE additional information is printed to the console.}}
    }
    \item{tikhonov:}{
      \itemize{
        \item{integration:}{'midpoint', 'trapezoidal', 'simpson'}
        \item{scale:}{'pdfstat', 'cdfstat', 'cdfpval'}
        \item{trim:}{0.01, 0.99}
        \item{symmetric:}{TRUE}
        \item{bin:}{'epdf', 'ecdf'}
        \item{from:}{-6}
        \item{to:}{6}
        \item{resolution:}{500}
        \item{method:}{'lcurve', 'gcv', 'aic'}
        \item{log:}{TRUE}
        \item{penalty:}{0}
        \item{lambda:}{10^seq(-10, 10, length=100)}
        \item{verbose:}{Default FALSE if TRUE additional information is printed to the console.}} }
    \item{'ferreira:'}{not yet implemeneted} }
}
\examples{
m <- 5000 ##number of genes
J <- 10 ##sample size per group
pi0 <- 0.8 ##proportion of non-differentially expressed genes
m0 <- as.integer(m*pi0)
mu <- rbitri(m - m0, a = log2(1.2), b = log2(4), m = log2(2)) #effect size distribution
data <- simdat(mu, m=m, pi0=pi0, J=J, noise=NULL)
library(genefilter)
stat <- rowttests(data, factor(rep(c(0, 1), each=J)), tstatOnly=TRUE)$statistic
pd <- pilotData(statistics=stat, samplesize=sqrt(J/2), distribution='norm')
ss <- sampleSize(pd, method='deconv')
plot(ss)
}
\author{
  Maarten van Iterson
}
\references{
  van Iterson, M., P. 't Hoen, P. Pedotti, G. Hooiveld, J. den Dunnen, G. van Ommen,
  J. Boer, and R. de Menezes (2009): 'Relative power and sample size analysis on
  gene expression profiling data,' BMC Genomics, 10, 439--449.

  Ferreira, J. and A. Zwinderman (2006a): 'Approximate Power and Sample Size
  Calculations with the Benjamini-Hochberg Method,' The International Journal
  of Biostatistics, 2, 1.

  Ferreira, J. and A. Zwinderman (2006b): 'Approximate Sample Size Calculations
  with Microarray Data: An Illustration,' Statistical Applications in Genetics and
  Molecular Biology, 5, 1.

  Hansen, P. (2010): Discrete Inverse Problems: Insight and Algorithms, SIAM: Fun-
  damentals of algorithms series.

  Langaas, M., B. Lindqvist, and E. Ferkingstad (2005): 'Estimating the proportion
  of true null hypotheses, with application to DNA microarray data,' Journal of the
  Royal Statistical Society Series B, 67, 555--572.

  Storey, J. (2003): 'The positive false discovery rate: A bayesian interpretation and
  the q-value,' Annals of Statistics, 31, 2013--2035.
}
\seealso{
  \code{\link{optim}}
}
