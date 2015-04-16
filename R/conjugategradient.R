##' Non-negative conjugate gradient algorithm
##'
##' C-implementation, details follow.
##' @title Non-negative conjugate gradient algorithm
##' @param A the A matrix of the system: Ax = b.
##' @param b the b vector of the system: Ax = b.
##' @param type for the conjugate-gradients method. Takes value '1' for the Fletcher-Reeves update, '2' for Polak-Ribiere and '3' for Beale-Sorenson.
##' @param trace tracing information on the progress of the optimization is produced.
##' @return list containg regression coefficients and some additional information.
##' @author Maarten van Iterson
nncg <- function(A, b, type=1, trace=FALSE)
  {

    if(type < 1 || type > 3)
      stop("Wrong type!")

    ##some error checking
    if(!all(is.finite(b)) || !all(is.finite(A)))
      stop("Data contains non-finite values!")

    if(ncol(A) != nrow(A) | nrow(A) != length(b))
      stop("Dimensions do not match!")

    output <- .C("nncg",
                 n = as.integer(ncol(A)),
                 x0 = as.vector(b*0, mode="double"),
                 A = as.vector(A, mode="double"),
                 b = as.vector(b, mode="double"),
                 fmin = as.double(0),
                 fail = as.integer(0),
                 type = as.integer(type),
                 trace = as.integer(trace),
                 objfcount = as.integer(0),
                 gradcount = as.integer(0), PACKAGE="SSPA")

    list(par=output$x0, value=output$fmin, counts=list('function'=output$objfcount, gradient=output$gradcount), convergence=output$fail)
  }

congrad <- function(object)
  {
    ##extract control parameters
    samplesize <- Samplesize(object)
    theta <- Theta(object)
    distribution <- distribution(object)
    bin <- object@control$bin
    scale <- object@control$scale
    trim <- object@control$trim
    symmetric <- object@control$symmetric
    verbose <- object@control$verbose
    from <- object@control$from
    to <- object@control$to

    if(scale == "cdfpval")
      statistics <- Pvalues(object)
    else
      statistics <- Statistics(object)

    resolution <- length(theta)+1

    tb <- trimbin(statistics, nbins=resolution, bin=bin, plot=verbose, symmetric=symmetric, trim=trim)

    x <- tb$x
    b <- tb$y

    ##extract functions from object
    df1 <- df1(object)
    pf1 <- pf1(object)
    qf0 <- qf0(object)

    N <- samplesize

    K <- switch(scale,
                pdfstat = function(x, y) df1(x=y, y=N*x),
                cdfstat = function(x, y) pf1(q=y, y=N*x),
                cdfpval = if(distribution == "norm" || distribution == "t")
                function(x, y) 1 - pf1(q=qf0(p=1 - y/2), y=N*x) + pf1(q=-qf0(p=1 - y/2), y=N*x)  ##two-sided t norm
                else
                function(x, y) 1 - pf1(q=qf0(p=1 - y), y=N*x) ##one-side f chisq
                )

    A <- midpoint(K, from, to, resolution-1, x)
    A <- cbind(K(0, x), t(A))

    beta <- nncg(A, b, trace=verbose)$par

    object@pi0 <- beta[1]
    if(beta[1] > 1)
      warning("Estimated pi0 > 1!")

    const <- sum(beta[-1]/(1 - beta[1])) * (theta[2] - theta[1])
    lambda.hat <- beta[-1]/(const * (1 - beta[1]))

    object@lambda <- lambda.hat
    object@theta <- theta
    
    object
  }
