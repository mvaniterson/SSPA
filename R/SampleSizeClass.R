setClass("SampleSize",
         representation("PilotData",
                        pi0           = "numeric",
                        lambda        = "numeric",
                        theta         = "numeric",
                        control = "list",
                        info    = "list"),
         prototype(PilotData     = "PilotData",
                   pi0           = numeric(1),
                   lambda        = numeric(1),
                   theta         = numeric(1),
                   control = list(),
                   info    = list())
         )

#######################################################################################
##Generic Accessor functions for class "SampleSize"
##
#######################################################################################
setGeneric("Pi0", function(object) { standardGeneric ("Pi0")})
setGeneric("Lambda", function(object) { standardGeneric ("Lambda")})
setGeneric("Theta", function(object) { standardGeneric ("Theta")})
setGeneric("Control", function(object) { standardGeneric ("Control")})
setGeneric("Info", function(object) { standardGeneric ("Info")})

setMethod("Pi0","SampleSize", function(object){ return(object@pi0)})
setMethod("Lambda","SampleSize", function(object){ return(object@lambda)})
setMethod("Theta","SampleSize", function(object){ return(object@theta)})
setMethod("Control","SampleSize", function(object){return(object@control)})
setMethod("Info","SampleSize", function(object){ return(object@info)})

#######################################################################################
##Show method for SampleSize
##
#######################################################################################
setMethod("show", signature("SampleSize"), function (object) { cat(str(object)) })

#######################################################################################
##Overload Generic plot method for class "SampleSize"
##plot density of effect sizes
#######################################################################################
setMethod("plot", signature(x="SampleSize"), definition = function(x, y, threshold = 0, ...)
          {
            object <- x

            dots <- list(...)
            main <- dots[["main"]]
            if(is.null(main)) main <- "density of effect sizes"
            xlab <- dots[["xlab"]]
            if (is.null(xlab)) xlab <- "effect size"
            ylab <- dots[["ylab"]]
            if (is.null(ylab)) ylab <- ""

            if(threshold != 0)
              lambda <- truncateEffectsize(object, threshold)
            else
              lambda <- Lambda(object)

            xyplot(lambda~Theta(object), type="l", xlab=xlab, ylab=ylab, main=main, ...)

          })

checking <- function(x, y)
{
  if(length(y) > 1)
    {
      if(is.character(y)[1])
        {
          if(length(setdiff(x, y)) > 0)
            stop(paste("Parameter '", x, "' is not allowed. Should be one of '", paste(y, collapse="', '"), "'!", sep=""))
        }
    }
  else
    {
      if(typeof(x) != typeof(y))
        stop(paste("type mismatch: ", x, y))
    }
}

deconvControl <- function(control)
{
  ##defaults
  defpar <- list(method = c("deconv","ferreira"),
                 pi0Method = c("Langaas", "Storey", "Ferreira", "Userdefined"),
                 pi0 = seq(0.1, 0.99, 0.01),
                 adjust = TRUE,
                 a=0.5,
                 bandwith = NULL,
                 kernel = c("fan", "wand", "sinc"),
                 from =-6,
                 to = 6,
                 resolution = 2^9,
                 verbose = FALSE)

  if(length(setdiff(names(control), names(defpar))) > 0)
    stop(paste("Unknown control parameter:", setdiff(names(control), names(defpar))), collapse=", ")

  ##set defaults
  conpar <- defpar

  ##overwrite defaults if given
  conpar[names(control)] <- control

  ##some more checks
  tmp <- lapply(names(conpar), function(x) checking(conpar[[x]], defpar[[x]]))

  ##get only one option
  conpar <- lapply(conpar, function(x){
    if(is.character(x[1]))
      x[1]
    else
      x
  })

  ##some logical checks
  if(conpar$pi0Method == "Userdefined" && all(names(conpar) != "pi0"))
    stop("must specify fraction of non-differently expressed genes!")
  else if(conpar$pi0Method == "Userdefined" && length(conpar$pi0) > 1)
    stop("must specify a single fraction of non-differently expressed genes!")
  else if(conpar$pi0Method == "Ferreira" && length(conpar$pi0) == 1)
    stop("must specify a grid of values!")
  else if(conpar$pi0 < 0 || conpar$pi0 > 1)
    stop("All pi0-values must be between 0 and 1")

  conpar
}

congradControl <- function(control)
{

  defpar <- list(method = "congrad",
                 integration = c("midpoint", "trapezoidal", "simpson"),
                 scale = c("pdfstat", "cdfstat", "cdfpval"),
                 trim = c(0.01, 0.99),
                 symmetric = TRUE,
                 bin = c("epdf", "ecdf"),
                 from = -6,
                 to = 6,
                 resolution = 500,
                 verbose = FALSE)

  if(length(setdiff(names(control), names(defpar))) > 0)
    stop(paste("Unknown control parameter:", setdiff(names(control), names(defpar))), collapse=", ")

  ##set defaults
  conpar <- defpar

  ##overwrite defaults if given
  conpar[names(control)] <- control

  ##some more checks
  tmp <- lapply(names(conpar), function(x) checking(conpar[[x]], defpar[[x]]))

  ##get only one option
  conpar <- lapply(conpar, function(x){
    if(is.character(x[1]))
      x[1]
    else
      x
  })

  ##some logical checks
  if(conpar$scale == "pdfstat")
    conpar$bin <- "epdf"
  else
    conpar$bin <- "ecdf"

  if(conpar$scale == "cdfpval")
    {
      conpar$symmetric <- FALSE
      conpar$trim <- c(0, 1)         #always between 0, 1 and sum to 1
    }

  if(conpar$trim[1] < 0 | conpar$trim[2] > 1)
    stop("trimming values should be within 0 and 1!")

  conpar
}

tikhonovControl <- function(control)
{
  ##defaults
  defpar <- list(method = "tikhonov",
                 integration = c("midpoint", "trapezoidal", "simpson"),
                 scale = c("pdfstat", "cdfstat", "cdfpval"),
                 trim = c(0.01, 0.99),
                 symmetric = TRUE,
                 bin = c("epdf", "ecdf"),
                 from = -6,
                 to = 6,
                 resolution = 500,
                 modelselection = c("lcurve", "gcv", "aic"),
                 log = TRUE,
                 penalty = 0,
                 lambda = 10^seq(-10, 10, length=100),
                 verbose = FALSE)

  if(length(setdiff(names(control), names(defpar))) > 0)
    stop(paste("Unknown control parameter:", setdiff(names(control), names(defpar))), collapse=", ")

  ##set defaults
  conpar <- defpar

  ##overwrite defaults if given
  conpar[names(control)] <- control

  ##some more checks
  tmp <- lapply(names(conpar), function(x) checking(conpar[[x]], defpar[[x]]))

  ##get only one option
  conpar <- lapply(conpar, function(x){
    if(is.character(x[1]))
      x[1]
    else
      x
  })

  ##check set binning
  if(conpar$scale == "pdfstat")
    conpar$bin <- "epdf"
  else
    conpar$bin <- "ecdf"


  conpar
}

defineEffectSizeRange <- function(object, from, to, resolution)
  {
    ##'TODO use trimmingbinning function for this!

    ##Resolution must be a power of 2, for the FFT and smaller than test statistics
    if(resolution >= length(Statistics(object)))
      {
        resolution <- 2^floor(log2(length(Statistics(object))))
        warning("Resolution should be smaller than number of test statistics!")
      }
    else if(resolution%%2 != 0)
      {
        resolution <- 2^ceiling(log2(resolution))
        warning("Resolution is set to a power of 2!")
      }

    ##Input range must be symmetric for the fast gnhat calculation
    if(abs(from) != to)
      {
        to <- max(c(abs(from), to))
        from = -to
        warning("Input range is made symmetric!")
      }

    seq(from = from, to = to, length = resolution)
  }

##' User friendly interface to class "SampleSize"
##'
##' @details
##' The default method is 'deconv' which is a kernel deconvolution density estimator implementated using fft.
##' The 'nncg' is a nonnegative conjugate gradient algorithm based on R's implementation see optim.
##' 'tikonov' implements ridge-regression with optimal penalty selection using the L-curve approach.
##' Higher order penalties are possible as well using a transformation to standard form.
##'  The 'control' argument is a list that can supply any of the following components, however per method some logical checks are built-in:
##' \itemize{
##'  \item{deconv:}{
##'  \itemize{
##'  \item{method:}{'deconv', 'ferreira'}
##'  \item{pi0Method:}{the pi0 estimation method one of 'Langaas', 'Storey', 'Ferreira', 'Userdefined'}
##'  \item{pi0:}{if method = 'ferreira' grid pi0-value need to be suppled e.g. seq(0.1, 0.99, 0.01)}
##'  \item{adjust:}{Default TRUE, adjust pi0 esitmate if density of effect size is somewhere negative.}
##'  \item{a:}{Adjust pi0 better approach suggested by Efron. Symmetric range around zero of size 0.5.}
##'  \item{bandwith:}{Default NULL uses  1/sqrt(log(length(statistics)))}
##'  \item{kernel:}{Either 'fan', 'wand', 'sinc' kernels can be used.}
##'  \item{from:}{Density of effect sizes should be estimated from = -6}
##'  \item{to:}{6}
##'  \item{resolution:}{Density of effect sizes should be estimated on 2^9 points.}
##'  \item{verbose:}{Default FALSE if TRUE additional information is printed to the console.}}
##' }
##'  \item{congrad:}{
##'  \itemize{
##'  \item{integration:}{'midpoint', 'trapezoidal', 'simpson'}
##'  \item{scale:}{'pdfstat', 'cdfstat', 'cdfpval'}
##'  \item{trim:}{0.01, 0.99}
##'  \item{symmetric:}{TRUE}
##'  \item{bin:}{'epdf', 'ecdf'}
##'  \item{from:}{-6}
##'  \item{to:}{6}
##'  \item{resolution:}{500}
##'  \item{verbose:}{Default FALSE if TRUE additional information is printed to the console.}}
##' }
##'  \item{tikhonov:}{
##'  \itemize{
##'  \item{integration:}{'midpoint', 'trapezoidal', 'simpson'}
##'  \item{scale:}{'pdfstat', 'cdfstat', 'cdfpval'}
##'  \item{trim:}{0.01, 0.99}
##'  \item{symmetric:}{TRUE}
##'  \item{bin:}{'epdf', 'ecdf'}
##'  \item{from:}{-6}
##'  \item{to:}{6}
##'  \item{resolution:}{500}
##'  \item{modelselection:}{'lcurve', 'gcv', 'aic'}
##'  \item{log:}{TRUE}
##'  \item{penalty:}{0}
##'  \item{lambda:}{10^seq(-10, 10, length=100)}
##'  \item{verbose:}{Default FALSE if TRUE additional information is printed to the console.}}
##' }
##'  \item{'ferreira:'}{}
##' }
##' @title User friendly interface to class 'SampleSize'
##' @param PilotData object of class 'PilotData'.
##' @param method estimation method one of 'deconv', 'congrad', 'tikhonov' or 'ferreira'. See 'Details'.
##' @param control A list of control parameters. See 'Details'.
##' @return object of class SampleSize.
##' @author Maarten van Iterson
##' @references
##' Langaas, Storey, Ferreira, Hansen, van Iterson
##' @seealso \code{\link{optim, fft}}
##' @export
##' @examples
##' m <- 5000 ##number of genes
##' J <- 10 ##sample size per group
##' pi0 <- 0.8 ##proportion of non-differentially expressed genes
##' m0 <- as.integer(m*pi0)
##' mu <- rbitri(m - m0, a = log2(1.2), b = log2(4), m = log2(2)) #effect size distribution
##' data <- simdat(mu, m=m, pi0=pi0, J=J, noise=NULL)
##' library(genefilter)
##' stat <- rowttests(data, factor(rep(c(0, 1), each=J)), tstatOnly=TRUE)$statistic
##' pd <- pilotData(statistics=stat, samplesize=sqrt(J/2), distribution='norm')
##' ss <- sampleSize(pd, method='deconv')
##' plot(ss)
sampleSize <- function(PilotData, method = c("deconv", "congrad", "tikhonov", "ferreira"), control=list(from=-6, to=6, resolution=2^9))
{
  ##create new SampleSize-object with PilotData
  object <- new("SampleSize", PilotData)

  ##verify method
  method <- match.arg(method)

  control[["method"]] <- method

  ##fill control parameters for the specific method
  control <- switch(method,
                    deconv = deconvControl(control),
                    congrad = congradControl(control),
                    tikhonov = tikhonovControl(control),
                    ferreira = deconvControl(control))

  if(method == "deconv")
    theta <- defineEffectSizeRange(object, control$from, control$to, control$resolution)
  else
    theta <- seq(from = control$from, to = control$to, length = control$resolution)

  ##update SampleSize-object
  object@control <- control
  object@theta <- theta

  ##estimate density of effect sizes
  object <- switch(method,
                   deconv = deconvolution(object),
                   congrad = congrad(object),
                   tikhonov = Tikhonov(object),
                   ferreira = Dn(object, doplot=control$verbose))
  object
}
