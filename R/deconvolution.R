##' TODO: suppress warning message in pt() :  full precision was not achieved in 'pnt'
##'TODO:replace massdist by trimbin?

##' deconvolution estimator using fft
##'
##' details follow
##' @title deconvolution estimator using fft
##' @param object of class SampleSize
##' @return object of class SampleSize
##' @author Maarten van Iterson
deconvolution <- function(object)
{

  ##extract control parameters
  statistics <- Statistics(object)
  theta <- Theta(object)
  samplesize <- Samplesize(object)

  ##get pi0 estimate or plugin value
  pi0 <- switch(object@control$pi0Method,
                Storey = qvalue(Pvalues(object))$pi0,
                Langaas = convest(Pvalues(object)),
                Userdefined = object@control$pi0)

  if(is.null(object@control$bandwidth))
    h <- 1/sqrt(log(length(statistics)))
  else
    h <- object@control$bandwidth

  adjust <- object@control$adjust
  a <- object@control$a
  kernel <- object@control$kernel

  ##start deconvolution of the mixture
  N <- 2*length(theta)

  from <- min(statistics) - 7 * h
  to <- max(statistics) + 7 * h

  ##make range symmetric around zero
  to <- max(c(abs(from), abs(to)))
  from <- sign(from)*to

  Dx <- 2*(to-from)/(N-1)
  Dt <- (2*pi)/(N*Dx)
  k <- n <- 0:(N-1)

  M <- N/2
  tn <- n*Dt
  tn[(M+2):N] <- -tn[M:2]

  xk <- k*Dx
  xk[(M+2):N] <- -xk[M:2]

  x <- seq.int(from, to, length.out = M)

  Nx <- length(statistics)

  ##empirical density estimation from library stat
  y <- .C("massdist",
          x = as.double(statistics),
          xmass = as.double(rep(1/Nx, Nx)),
          nx = as.integer(Nx),
          xlo = as.double(from),
          xhi = as.double(to),
          y = double(N),
          ny = as.integer(M),
          PACKAGE = "SSPA")$y

  f0n <- c(df0(object)(x), rep(0, M))*Dx
  cf0n <- fft(f0n) ##  cf0(object)(tn)
  
  ##check if density goes below zero adjust pi0
  diff <- y - pi0*f0n
  if(any(diff < -1e-5) & adjust)        ##default -1e-5
    {
      object@info$pi0 <- pi0 ##store unadjusted estimated pi0
      if(a == 0)
        {
          idx <- which.min(diff)
          pi0 <- y[idx]/f0n[idx]
        }
      else
        {
          ##Better Upper Bound Estimates for pi0 Efron et al. (2001) Empirical Bayes Analysis of a Microarray Experiment
          indices <- which(x > -a & x < a)
          pi0 <- sum(y[indices])/sum(f0n[indices]) ##or 2*pnorm(a) - 1 or 2 *pt(a, df=obj@dof)- 1
        }
    }

  fn <- (y - pi0*f0n)/(1 - pi0)
  
  epsilon <- 1e-10 ##circumvent division by zero
  lambda <- fft(fft(fn)*Conj(Kernel(tn, h, kernel)/Dx)/(cf0n+epsilon), inverse=TRUE)

  lambda <- pmax.int(0, Re(lambda))/N

  ##scale since int k(x, ay)g(y)dy by change of variable z =  ay gives int k(x, z) g(z/a)/a dz
  ##lambda on user-defined scale
  lambda <- approx(xk/samplesize, lambda*samplesize, theta, rule=2)$y

  ##normalization
  const <- sum(lambda)*(theta[2]-theta[1])

  object@lambda <- lambda/const
  object@theta <- theta
  
  object@pi0 <- pi0 ##store adjusted pi0

  object
}

#######################################################################################
##truncate density of effect-sizes for given threshold
#######################################################################################
truncateEffectsize <- function(object, threshold)
{
  index.upper <- which.min(abs(Theta(object)-threshold))
  lambda <- Lambda(object)
  lambda[(length(Theta(object))-index.upper + 1):index.upper] <- 0.0
  lambda/(sum(lambda)*(Theta(object)[2]-Theta(object)[1]))
}
