#######################################################################################
##Class SampleSize
#######################################################################################
setClass("SampleSize",
	representation("PilotData",
                       pi0       = "list",								 
                       lambda    = "numeric",
                       theta     = "numeric", 
                       adjust    = "logical",		
                       method    = "character",
                       bandwidth = "numeric",
                       kernel    = "character",
                       nKnots 	 = "numeric", 
                       bDegree   = "numeric"),
	prototype(PilotData = "PilotData",                  
            pi0       = list(),
            lambda    = numeric(1),
            theta     = numeric(1),
            adjust    = logical(1),		
            method    = character(1),
            bandwidth	= numeric(1),
            kernel    = character(1),
            nKnots    = numeric(1), 
            bDegree   = numeric(1))
)

#######################################################################################
##Generic Accessor functions for class "SampleSize"
#
#######################################################################################
setGeneric("pi0", function(object) { standardGeneric ("pi0")})
setGeneric("lambda", function(object) { standardGeneric ("lambda")})
setGeneric("theta", function(object) { standardGeneric ("theta")})

setMethod("pi0","SampleSize", function(object){ return(object@pi0)})
setMethod("lambda","SampleSize", function(object){ return(object@lambda)})
setMethod("theta","SampleSize", function(object){ return(object@theta)})

#######################################################################################
##Show method for SampleSize
#
#######################################################################################
setMethod("show", signature("SampleSize"), function (object)
{

  cat("An object of class \"", class(object),"\"\n", sep = "")	
  cat("Distribution of effect sizes is estimated from ", object@theta[1], " to ", 
      object@theta[length(object@theta)], " using ", length(object@theta), 
      " points.", "\n", sep="")
  cat("Method for estimation of the proportion of non-differentially expressed: \"", object@method, "\n", sep = "")
  if(object@method!="Ruppert")	
    {
      Gm <- ifelse(length(object@pi0) == 2, object@pi0[[2]], 0.0)		
      cat("Fraction of non-differentially expressed genes: ", signif(object@pi0[[1]], 4), 
          " (adjusted=", signif(Gm, 4), ").\n", sep = "")
      cat("Kernel used in the deconvolution is \"", object@kernel, "\" with bandwidth ", 
          signif(object@bandwidth, 4), ".\n", sep = "")
    }	else {
      cat("Fraction of non-differentially expressed genes: \n", 
          names(object@pi0)[1], " = ", signif(object@pi0[[1]], 4), ".\n", sep = "")
      cat(names(object@pi0)[2], " = ", signif(object@pi0[[2]], 4), ".\n", sep = "")
      cat(names(object@pi0)[3], " = ", signif(object@pi0[[3]], 4), ".\n", sep = "")
      cat("Bsplines of degree ", object@bDegree, " with ", object@nKnots, " knots are used.\n", sep = "")
    }     

})


#######################################################################################
##Overload Generic plot method for class "SampleSize" 
#plot density of effect sizes
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
    lambda <- lambda(object)

  xyplot(lambda~theta(object), type="l", xlab=xlab, ylab=ylab, main=main, ...)		

})

#######################################################################################
##User friendly interface to SampleSize
#
#######################################################################################
sampleSize <- function(PilotData, method = c("Langaas", "Storey", "Ferreira", "Ruppert", "Userdefined"), from = -6, to = 6, resolution = 2^10, kernel = c("fan", "wand", "sinc"), pi0 = seq(0.1, 0.99, 0.01), adjust = TRUE, a=0.5, bandwidth, nKnots = 11, bDegree = 3, ...) 
{
  ##create new SampleSize-object with PilotData
  object <- new("SampleSize", PilotData)
	
  ##Resolution must be a power of 2, for the FFT and smaller than test statistics
  if(resolution >= length(statistics(PilotData)))
    {
      resolution <- 2^floor(log2(length(statistics(PilotData))))
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

  object@theta <- seq(from = from, to = to, length = resolution)	

  ##estimate fraction of p-values computed under null hypothesis
  object@method <- match.arg(method)	

  if(object@method == "Userdefined" && missing(pi0))  
    stop("must specify fraction of non-differently expressed genes!")
  else if(object@method == "Userdefined" && length(pi0) > 1) 
    stop("must specify a single fraction of non-differently expressed genes!")
  else if(object@method == "Ferreira" && length(pi0) == 1)
    stop("must specify a grid of values!")
  else if(any(pi0 < 0 || pi0 > 1))     #actually 0 and 1 make no sense
    stop("All pi0-values must be between 0 and 1")

  object@adjust <- adjust

  if(object@method!="Ruppert") {
    if(missing(bandwidth))
      object@bandwidth <- 1/sqrt(log(length(statistics(object)))) #default bandwidth
    else
      object@bandwidth <- bandwidth
    object@kernel <- match.arg(kernel)	
    object@pi0[[object@method]] <- switch(object@method, Storey = qvalue(pvalues(object), ...)$pi0,
                                          Langaas = convest(pvalues(object), ...),
                                          Ferreira = Dn(object, adjust, a, pi0, ...),		
                                          Userdefined = pi0)		
		
    object <- Deconvolution(object, adjust, a)		
		
  } else {			
    object@nKnots <- nKnots
    object@bDegree <- bDegree
    object <- psreg(object)		
  }
  object
}  

##estimate density of effect sizes
#new
deconvolution <- function(object, adjust, a)
{
		
  N <- 2*length(theta(object))	
  h <- object@bandwidth

  from <- min(statistics(object)) - 7 * h
  to <- max(statistics(object)) + 7 * h

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

  Nx <- length(statistics(object))

  ##empirical density estimation from library stat
  y <- .C("massdist", 
          x = as.double(statistics(object)), 
          xmass = as.double(rep(1/Nx, Nx)),
          nx = as.integer(Nx), 
          xlo = as.double(from), 
          xhi = as.double(to), 
          y = double(N),
          ny = as.integer(M), 
          PACKAGE = "SSPA")$y
		
  cnorm <- function(t, mean, sd) exp(complex(imaginary = mean*t, real = - 0.5*(sd*t)^2))

  ##from QRMlib or prob package:
  ct <- function(t, df){
    aux <- ifelse(t == 0, 0, log(besselK(abs(t)*sqrt(df), df/2, expon.scaled=FALSE) * ((abs(t)*sqrt(df))^(df/2))*(2^(1-df/2))) -lgamma(df/2))
    return(exp(aux))
  }

  if(distribution(object) == "Normal")
    {	
      null <- function(x) dnorm(x, 0, 1) 		
      cnull <- function(x) cnorm(x, 0, 1)
    } else if(distribution(object) == "Student"){
      null <- function(x) dt(x, df=dof(object))
      cnull <- function(x) ct(x, df=dof(object))
    } 
	
  f0n <- c(null(x), rep(0, M))*Dx
  cf0n <- fft(f0n) ##cnull(tn)		
	
  ##check if density goes below zero adjust pi0	 
  Gm <- object@pi0[[1]]	
  diff <- y - Gm*f0n
  if(any(diff < -1e-5) & adjust)        #default -1e-5
    {	
      if(a == 0) 
        {
          idx <- which.min(diff)
          Gm <- y[idx]/f0n[idx]
          object@pi0[["Adjusted"]] <- Gm		
        }
      else 
        {	
                                        ##Better Upper Bound Estimates for pi0 Efron et al. (2001) Empirical Bayes Analysis of a Microarray Experiment		
          indices <- which(x > -a & x < a)
          Gm <- sum(y[indices])/sum(f0n[indices]) #or 2*pnorm(a) - 1 or 2 *pt(a, df=obj@dof)- 1
          object@pi0[["Adjusted"]] <- Gm
        }		
    }	
		
  kernel <- switch(object@kernel, fan = ifelse(abs(h*tn) > 1, 0, (1-(h*tn)^2)^3),
                   wand = ifelse(abs(h*tn) < 1/2 | abs(h*tn) > 1, 0, (2*(1-abs(h*tn))^3)) + 
                   ifelse(abs(h*tn) > 1/2, 0, (1-6*(h*tn)^2+6*abs(h*tn)^3)),
                   sinc = ifelse(abs(h*tn) > 1, 0, (1-abs(h*tn))))/Dx
		
  fn <- (y - Gm*f0n)/(1 - Gm) 
	
  epsilon <- 1e-10                      #circumvent division by zero
  lambda <- fft(fft(fn)*Conj(kernel)/(cf0n+epsilon), inverse=TRUE)*sqrt(effectivesamplesize(object))		
		
  lambda <- pmax.int(0, Re(lambda))/N 
	
  lambda <- approx(xk/sqrt(effectivesamplesize(object)), lambda, theta(object), rule=2)$y	 

                                        ##normalization
  const <- sum(lambda)*(theta(object)[2]-theta(object)[1])
  object@lambda <- lambda/const	

  object
}

##estimate density of effect sizes
#old
##estimate density of effect sizes
Deconvolution <- function(object, adjust, a)
{
  N <- length(theta(object))
  M <- length(statistics(object))  	
  h <- object@bandwidth
	
  ##Ferreira suggested 
  ##4*max(abs(obj@testStatistics))
	
  lo <- (min(statistics(object)) - 7*h) #default 7
  up <- (max(statistics(object)) + 7*h)
	
  delta <- 2*(up-lo)/(2*N-1)
	
  tk <- (1:(2*N))*delta
  tk[(N+2):(2*N)] <- -tk[N:2]

  xk <- (1:(2*N ))*(2*pi)/(delta*2*N)
  xk[(N+2):(2*N)] <- -xk[N:2]

  x <- seq(from = lo, to = up, length = N)

  ##empirical density estimation from library stat
  y <- .C("massdist", 
          x = as.double(statistics(object)), 
          xmass = as.double(rep(1/M, M)),
          nx = as.integer(M), 
          xlo = as.double(lo), 
          xhi = as.double(up), 
          y = double(2*N),
          ny = as.integer(N), 
          PACKAGE = "SSPA")$y
	
  ##density of effect-sizes	
  fz <- c(delta*dnorm(x, 0, 1), numeric(N)) 

  if(distribution(object) == "Normal")
    {	
      ##characteristic function of Normal	
      Fz <- fft(dnorm(tk, 0, 1))		
    } else {
      ##characteristic function of Student t
      Fz <- fft(dt(tk, df=dof(object)))
    }
	
  ##check if density goes below zero adjust pi0	 
  Gm <- object@pi0[[1]]	
  diff <- y - Gm*fz
  if(any(diff < -1e-5) & adjust)        #default -1e-5
    {	
      idx <- which.min(diff)
      Gm <- y[idx]/fz[idx]
      object@pi0[["Adjusted"]] <- Gm		
    }
	
  y <- (y - Gm*fz)/(1-Gm)				
	
  ##add reference to kernels
  Kz <- switch(object@kernel, fan = ifelse(abs(h*xk) > 1, 0, (1-(h*xk)^2)^3),
               wand = ifelse(abs(h*xk) < 1/2 | abs(h*xk) > 1, 0, (2*(1-abs(h*xk))^3)) + 
               ifelse(abs(h*xk) > 1/2, 0, (1-6*(h*xk)^2+6*abs(h*xk)^3)),
               sinc = ifelse(abs(h*xk) > 1, 0, (1-abs(h*xk))))
	
  V <- Kz/(delta*Conj(Fz))
  V[is.na(V)] <- 0.0                    #maybe nicer solution

  ##deconvolution
  lambda <- fft(fft(y)*V, inverse = TRUE)*sqrt(effectivesamplesize(object))
	
  ##normalize
  lambda <- Re(lambda)[1:N]/(2*N)
	
  ##user defined scale
  lambda <- approx(x/sqrt(effectivesamplesize(object)), lambda, theta(object), rule = 2)$y #rule=2 no NA's are introduced
	
  ##normalization
  const <- sum(lambda)*(theta(object)[2]-theta(object)[1])
  object@lambda <- lambda/const	

  object
}




##parametric estimation of the distribution of p-values computed under the alternative hypotheses
Gnhat <- function(object, effectivesamplesize, threshold = 0, lower = 0, upper = 1.0, resolution)
{	
  u <- seq(from = lower + .Machine$double.eps, to = upper - .Machine$double.eps, length = resolution) 

  if(missing(effectivesamplesize)) effectivesamplesize <- effectivesamplesize(object)

  if(threshold != 0)
    lambda <- truncateEffectsize(object, threshold)
  else
    lambda <- lambda(object)

  theta <- theta(object)
  ##theta is symmetric thus time can be reduced by a half	
  ##theta <- seq(0, max(theta(object)), length=length(theta(object))/2)
		
  if(distribution(object)=="Normal")
    {	
      ##normal case	
      gnhat <- outer(qnorm(1-u/2, 0, 1), theta*sqrt(effectivesamplesize), function(x, y) 1 - pnorm(x-y) + pnorm(-x-y))				
    } else if(distribution(object)=="Student"){
      ##student case THIS NOT "the shifted-t" but the exact noncentral t
      ##TODO: suppress warning message in pt(-qt(1 - x/2, df = obj@dof), df = obj@dof, ncp = y) ... :
                                        #                                  full precision was not achieved in 'pnt'
      ## and using the t-dist is slow	
      gnhat <- outer(qt(1-u/2, df=dof(object), ncp=0), theta*sqrt(effectivesamplesize),
                     function(x, y) 1 - pt(x, df=dof(object), ncp=y) + pt(-x, df=dof(object), ncp=y))
    }
	
  ##theta is symmetric thus time can be reduced by a half	
  ##gnhat <- cbind(gnhat[,length(theta):1], gnhat[,1:length(theta)])
	
  gnhat%*%(lambda*(theta[2]-theta[1]))
}

##Bernstein estimate of the distribution function of p-values
Hntilde <- function(pvalues, resolution)
{
  n <- length(pvalues)
  x <- seq(from = 0, to = 1, length=resolution)
  m <- round(n^(2/3))	
  edf <- sapply(1:m, function(x) length(pvalues[pvalues <= x/m])/n, USE.NAMES = FALSE)		
  return(rowSums(matrix(edf*dbinom(1:m, m, rep(x, each = m)), ncol = m, byrow = TRUE)))
}

##non-parametric estimate of the distribution of p-values computed under the alternative hypotheses
Gntilde <- function(H.n.tilde, pi0)
{
  F <- seq(from = 0, to = 1, length = length(H.n.tilde)) #uniform distribution
  return((H.n.tilde - F*pi0)/(1 - pi0))
}

##parametric estimation of fraction on non-differently expressed genes( according to Ferreira)
Dn <- function(object, adjust, a, pi0, xlim = c(0,1), resolution = 1000, doplot = FALSE)
{	
  H.n.tilde <- Hntilde(pvalues(object), resolution)
  Dn <- numeric(length(pi0))
  for(i in seq(along = pi0))
    {
      object@pi0[["Ferreira"]] <- pi0[i]
      x <- deconvolution(x, adjust, a)	
      G.n.tilde <- Gntilde(H.n.tilde, pi0[i])
      G.n.hat <- Gnhat(x, lower = xlim[1], upper = xlim[2], resolution = resolution)		
      Dn[i] <- sum(abs(G.n.hat - G.n.tilde))/length(G.n.tilde) #distance measure
    }	
  x <- which.min(Dn)
  y <- pi0[x]

  if(doplot)
    {
      plot(pi0, Dn, type = 'b')      #maybe store in samplesize object
      points(x, y, pch = 19)	
    }
  y
}

##truncate density of effect-sizes for given threshold
truncateEffectsize <- function(object, threshold)
{
  index.upper <- which.min(abs(theta(object)-threshold))
  lambda <- lambda(object)
  lambda[(length(theta(object))-index.upper + 1):index.upper] <- 0.0
  lambda/(sum(lambda)*(theta(object)[2]-theta(object)[1]))	
}

##power estimation
power <- function(object, threshold = 0, fdr = 0.1, effectivesamplesizes = NULL, plot=FALSE, type = 'l', ylim = c(0,1), xlim = c(0,1), 
                  xlab = "p-value", ylab = "average power", main="Power curve(s)", sub, ...)
{	
  resolution <- 1000                    #default
  u <- as.matrix(seq(from = xlim[1], to = xlim[2], length = resolution))	

  if(object@method=="Ruppert")
    Gm <- object@pi0[[2]]               # use semi compromise
  else
    Gm <- ifelse(length(object@pi0) == 1, object@pi0[[1]], object@pi0[[2]]) #use adjusted pi0 if calculated

                                        #Note by Ferreira 
  if(all(fdr >= Gm))
    warning("fdr should be smaller than the proportion of non-differentially expressed genes!")

  r.hat <- fdr*(1 - Gm)/(Gm*(1 - fdr))

  if(is.null(effectivesamplesizes))
    effectivesamplesizes <- effectivesamplesize(object) 
	
  G.n.hat <- matrix(numeric(), nrow = resolution, ncol = length(effectivesamplesizes))		
  for(i in seq(along = effectivesamplesizes))
    G.n.hat[,i] <- Gnhat(object, effectivesamplesizes[i], threshold=threshold, lower = xlim[1], upper = xlim[2], resolution = resolution)
 
  ##maybe different way to estimate the intersect
  power <- uStar <- matrix(numeric(length(effectivesamplesizes)*length(fdr)), nrow=length(effectivesamplesizes))
  for(i in 1:length(effectivesamplesizes))
    {
      for(j in 1:length(fdr))
        {
          id <- which(G.n.hat[,i] > u/r.hat[j])	
          if(length(id) <= 1)
            {
              uStar[i,j] <-	0.0001
              power[i,j] <- 0
            }
          else
            {
              id <- max(id)
              uStar[i,j] <-	0.5*(u[id] + u[id-1])
              power[i,j] <- 0.5*(G.n.hat[id,i] + G.n.hat[id-1,i])			
            }
        }
    }
	
  ##plot the power
  if(plot)
    {
      matplot(u, G.n.hat, type = type, xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab, 
              col=1, lty=1, ...)					    

      pch.power <- 1:6 
      lty.fdr <- 1:6
      pch.legend <- c(rep(NA, length(fdr)), as.numeric(pch.power))
      lty.legend <- c(lty.fdr[1:length(fdr)], 1)

      pts <- seq(1, resolution, 25)
		
      matlines(u[pts], G.n.hat[pts,], type = 'p', col=1, pch=pch.power)
      matlines(u, tcrossprod(u, 1/as.matrix(r.hat)), lty=lty.fdr, col=1)

      x0 <- numeric(nrow(power)*ncol(power))
      x1 <- as.vector(uStar)
      y0 <- y1 <- as.vector(power)

      arrows(x0, y0, x1, y1, code = 1, length = 0.05, col=1, lty=rep(lty.fdr, each=ncol(G.n.hat))) #not really flexible
      text(x1 + xlim[2]/12.5, y1, signif(power, 3), col=1)		
      legend("bottomright", c(paste("fdr", fdr), paste("sample size", round(2*effectivesamplesizes, 2))), bty='n', lty=lty.legend, pch=pch.legend)

    }	
  colnames(power) <- paste("fdr =", fdr)
  rownames(power) <- paste("sample size =", round(2*effectivesamplesizes, 2)) #remove the two

  power
}


