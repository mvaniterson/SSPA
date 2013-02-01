##' bisection method for finding ustar
##'
##' solve G(u) = u pi0 (1-alpha)/(aplha(1-pi0)), more  details follow.
##' @title bisection method for finding ustar
##' @param g left-hand-side
##' @param h right-hand side
##' @param umax between [0, 1] usually alpha.
##' @param verbose TRUE/FALSE
##' @param plot TRUE/FALSE
##' @return ustar
##' @author Maarten van Iterson
findroot <- function(g, h, umax, verbose=FALSE, plot=FALSE)
  {
    if(plot)
      {
        curve(g, 0, 1, n=1000, lwd=2, xlab="u", ylab=expression(hat(G)[1](u)), xlim=c(0, 0.1))
        curve(h, add=TRUE, lwd=2, lty=2)
        grid()
      }

    u <- umax

    while(abs(g(u) - h(u)) > 0.01)
      {
        u <- u + sign(g(u) - h(u))*u/2
        if(u < 0 | u > 1)
          {
            print("u outside boundary!!!") ##or return -Inf or +Inf
            break;
          }
        if(verbose)
          print(u)

        if(plot)
          {
            abline(v=u, col=1, lty=3)
            abline(h=g(u), col=1, lty=3)
          }
      }

    if(plot)
      {
        curve(g, add=TRUE, lwd=2, n=1000)
        curve(h, add=TRUE, lwd=2, lty=2)
        abline(v=u, col=2)
        abline(h=g(u), col=2)
      }
    u
  }

##' Predict power for given vector of sample sizes
##'
##' details follow.
##' @title Predict power for given vector of sample sizes
##' @param object of class 'SampleSize'
##' @param samplesizes vector of total sample sizes.
##' @param alpha FDR.
##' @param verbose TRUE/FALSE
##' @param plot TRUE/FALSE
##' @return predicted power.
##' @author Maarten van Iterson
##' @export
predictpower <- function(object, samplesizes, alpha=0.1, verbose=FALSE, plot=FALSE)
{
  lambda <- Lambda(object)
  theta <- Theta(object)
  pi0 <- Pi0(object)
  distribution <- distribution(object)

  if(!is.numeric(alpha))
    stop("FDR-level shoud be numeric!")
  if(alpha <= 0 || alpha >= 1)
    stop("FDR-level should be between [0,1]!")  
  if(alpha >= pi0)   ##note by Ferreira
    warning("FDR-level should be smaller than the proportion of non-differentially expressed genes!")
  
  ##extract functions from object
  pf1 <- pf1(object)
  qf0 <- qf0(object)

  if(distribution == "norm" || distribution == "t")
    G <- function(x, y, N) 1 - pf1(q=qf0(p=1 - y/2), y=N*x) + pf1(q=-qf0(p=1 - y/2), y=N*x)  ##two-sided t norm
  else
    G <- function(x, y, N) 1 - pf1(q=qf0(p=1 - y), y=N*x) ##one-side f chisq

  if(missing(samplesizes))
    samplesizes <- Samplesize(object)

  h <- function(x, a=pi0*(1-alpha)/(alpha*(1-pi0))) a*x
  w <- lambda*(theta[2]-theta[1])
  averagePower <- numeric(length(samplesizes))
  for(i in 1:length(samplesizes))
    {
      ##midpoint
      g <- function(u, t=theta, n=samplesizes[i], w = lambda*(theta[2]-theta[1])) sapply(u, function(x) sum(w*G(t, x, n)))
      u <- findroot(g, h, alpha, verbose=verbose, plot=plot)
      averagePower[i] <- g(u)
    }
  averagePower
}









