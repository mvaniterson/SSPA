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

##Show method for SampleSize
setMethod("show", signature("SampleSize"), 
function (object) {

	cat("An object of class \"", class(object),"\"\n", sep = "")	
	cat("Distribution of effect sizes is estimated from ", object@theta[1], " to ", 
      object@theta[length  (object@theta)], " using ", length(object@theta), 
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

##User friendly interface to SampleSize
sampleSize <- function(PilotData, method = c("Langaas", "Storey", "Ferreira", "Ruppert", "Userdefined"), from = -6, to = 6, resolution = 2^10, kernel = c("fan", "wand", "sinc"), pi0 = seq(0.1, 0.99, 0.01), adjust = TRUE, nKnots = 11, bDegree = 3, ...) 
{
	#create new SampleSize-object with PilotData
	obj <- new("SampleSize", PilotData)
	
	#Resolution must be a power of 2, for the FFT and smaller than test statistics
	if(resolution >= length(PilotData@testStatistics))
	{
		resolution <- 2^floor(log2(length(PilotData@testStatistics)))
    warning("Resolution should be smaller than number of test statistics!")    
	}
  else if(resolution%%2 != 0)
  {
		resolution <- 2^ceiling(log2(resolution))
    warning("Resolution is set to a power of 2!")
	}
	
	#Input range must be symmetric for the fast gnhat calculation
	if(abs(from) != to)
	{
		to <- max(c(abs(from), to))
		from = -to
		warning("Input range is made symmetric!")
	}

	obj@theta <- seq(from = from, to = to, length = resolution)	

	#estimate fraction of p-values computed under null hypothesis
	obj@method <- match.arg(method)	

	if(obj@method == "Userdefined" && missing(pi0))  
  	stop("must specify fraction of non-differently expressed genes!")
	else if(obj@method == "Userdefined" && length(pi0) > 1) 
  	stop("must specify a single fraction of non-differently expressed genes!")
	else if(obj@method == "Ferreira" && length(pi0) == 1)
  	stop("must specify a grid of values!")
	else if(any(pi0 < 0 || pi0 > 1)) #actually 0 and 1 make no sense
    stop("All pi0-values must be between 0 and 1")

	obj@adjust <- adjust

	if(obj@method!="Ruppert") {
		obj@bandwidth <- 1/sqrt(log(length(obj@testStatistics))) #default bandwidth
		obj@kernel <- match.arg(kernel)	
		obj@pi0[[obj@method]] <- switch(obj@method, Storey = qvalue(obj@pValues, ...)$pi0,
                                                Langaas = convest(obj@pValues, ...),
   																				      Ferreira = Dn(obj, adjust, pi0, ...),
			  	                                      Userdefined = pi0)		
		obj <- deconvolution(obj, adjust)	
	} else {			
		obj@nKnots <- nKnots
		obj@bDegree <- bDegree
		obj <- psreg(obj)		
	}
	obj
}  

##estimate density of effect sizes
deconvolution <- function(obj, adjust)
{
	N <- length(obj@theta)
	M <- length(obj@testStatistics)  	
	h <- obj@bandwidth
	
	#Ferreira suggested 
	#4*max(abs(obj@testStatistics))
	
  lo <- (min(obj@testStatistics) - 7*h) #default 7
	up <- (max(obj@testStatistics) + 7*h)
	
	delta <- 2*(up-lo)/(2*N-1)
	
	tk <- (1:(2*N))*delta
	tk[(N+2):(2*N)] <- -tk[N:2]

	xk <- (1:(2*N ))*(2*pi)/(delta*2*N)
	xk[(N+2):(2*N)] <- -xk[N:2]

	x <- seq(from = lo, to = up, length = N)

	#empirical density estimation from library stat
	y <- .C("massdist", 
          x = as.double(obj@testStatistics), 
          xmass = as.double(rep(1/M, M)),
          nx = as.integer(M), 
          xlo = as.double(lo), 
          xhi = as.double(up), 
          y = double(2*N),
				  ny = as.integer(N), 
          PACKAGE = "stats")$y
	
	
	#TODO: Maybe replace this with analytical characteristic functions
	cnorm <- function(t, mean, sd) exp(complex(imaginary = mean*t, real = - 0.5*(sd*t)^2))
	# based on QRMlib package:
	ct <- function(t, df){
    aux <- ifelse(t == 0, 0, log(besselK(abs(t)*sqrt(df), df/2, expon.scaled=FALSE) * ((abs(t)*sqrt(df))^(df/2))*(2^(1-df/2))) -lgamma(df/2))
		return(exp(aux))
	}



	if(obj@nullDist == "normal")
	{	
		#density of effect-sizes	
		fz <- c(delta*dnorm(x, 0, 1), numeric(N)) 
		#characteristic function of Normal	
		#Fz <- fft(dnorm(tk, 0, 1))
		Fz <- cnorm(xk, 0, 1)			
	} else {
		fz <- c(delta*dt(x, df=obj@dof), numeric(N)) 
		#characteristic function of Student t
		#Fz <- fft(dt(tk, df=obj@dof))
		Fz <- ct(xk, df=obj@dof)			
	} 
	
	#check if density goes below zero adjust pi0	 
	Gm <- obj@pi0[[1]]	
	diff <- y - Gm*fz
	if(any(diff < -1e-5) & adjust) #default -1e-5
	{	
		#idx <- which.min(diff)
		#Gm <- y[idx]/fz[idx]
		#obj@pi0[["Adjusted"]] <- Gm		
		#Better Upper Bound Estimates for pi0 Efron et al. (2001) Empirical Bayes Analysis of a Microarray Experiment
		a <- 0.5	
	  indices <- which(x > -a & x < a)
		Gm <- sum(y[indices])/sum(fz[indices]) #or 2*pnorm(a) - 1 or 2 *pt(a, df=obj@dof)- 1
		obj@pi0[["Adjusted"]] <- Gm		
	}
	
	y <- (y - Gm*fz)/(1-Gm)				
	
	#TODO: add references to kernels
	Kz <- switch(obj@kernel, fan = ifelse(abs(h*xk) > 1, 0, (1-(h*xk)^2)^3),
                           wand = ifelse(abs(h*xk) < 1/2 | abs(h*xk) > 1, 0, (2*(1-abs(h*xk))^3)) + 
                                  ifelse(abs(h*xk) > 1/2, 0, (1-6*(h*xk)^2+6*abs(h*xk)^3)),
                           sinc = ifelse(abs(h*xk) > 1, 0, (1-abs(h*xk))))
	
  V <- Kz/(delta*Conj(Fz))
	V[is.na(V)] <- 0.0 #maybe nicer solution

	#deconvolution
	lambda <- fft(fft(y)*V, inverse = TRUE)*sqrt(obj@sampleSize)
	
	#normalize
  lambda <- Re(lambda)[1:N]/(2*N)
	
	#user defined scale
	lambda <- approx(x/sqrt(obj@sampleSize), lambda, obj@theta, rule = 2)$y #rule=2 no NA's are introduced
	
  #normalization
	const <- sum(lambda)*(obj@theta[2]-obj@theta[1])
	obj@lambda <- lambda/const	

	obj
}

##parametric estimation of the distribution of p-values computed under the alternative hypotheses
Gnhat <- function(obj, sampleSize, threshold = 0, lower = 0, upper = 1.0, resolution)
{	
	u <- seq(from = lower + .Machine$double.eps, to = upper - .Machine$double.eps, length = resolution) 

	if(missing(sampleSize)) sampleSize <- obj@sampleSize

	if(threshold!=0)
		lambda <- truncateEffectsize(obj, threshold)
	else
		lambda <- obj@lambda

	if(obj@nullDist=="normal")
	{	
		#normal case	
		gnhat <- outer(qnorm(1-u/2, 0, 1), obj@theta*sqrt(sampleSize), function(x, y) 1 - pnorm(x-y) + pnorm(-x-y))				
	} else {
		#student case
		#TODO: suppress warning message in pt(-qt(1 - x/2, df = obj@dof), df = obj@dof, ncp = y) ... :
    #                                  full precision was not achieved in 'pnt'
		# and using the t-dist is slow	
		# gnhat <- outer(qt(1-u/2, df=obj@dof, ncp=0), obj@theta*sqrt(sampleSize),
    #                 function(x, y) 1 - pt(x, df=obj@dof, ncp=y) + pt(-x, df=obj@dof, ncp=-y))
		gnhat <- outer(qt(1-u/2, df=obj@dof, ncp=0), obj@theta*sqrt(sampleSize), function(x, y) 2*(1 - pt(x, df=obj@dof, ncp=y)))
	}


	gnhat%*%(lambda*(obj@theta[2]-obj@theta[1]))
}

##Bernstein estimate of the distribution function of p-values
Hntilde <- function(pValues, resolution)
{
	n <- length(pValues)
	x <- seq(from = 0, to = 1, length=resolution)
	m <- round(n^(2/3))	
	edf <- sapply(1:m, function(x) length(pValues[pValues <= x/m])/n, USE.NAMES = FALSE)		
	return(rowSums(matrix(edf*dbinom(1:m, m, rep(x, each = m)), ncol = m, byrow = TRUE)))
}

##non-parametric estimate of the distribution of p-values computed under the alternative hypotheses
Gntilde <- function(H.n.tilde, pi0)
{
	F <- seq(from = 0, to = 1, length = length(H.n.tilde)) #uniform distribution
	return((H.n.tilde - F*pi0)/(1 - pi0))
}

##parametric estimation of fraction on non-differently expressed genes( according to Ferreira)
Dn <- function(x, adjust, pi0, xlim = c(0,1), resolution = 1000, doplot = FALSE)
{	
	H.n.tilde <- Hntilde(x@pValues, resolution)
	Dn <- numeric(length(pi0))
	for(i in seq(along = pi0))
	{
		x@pi0[["Ferreira"]] <- pi0[i]
		x <- deconvolution(x, adjust)	
		G.n.tilde <- Gntilde(H.n.tilde, pi0[i])
		G.n.hat <- Gnhat(x, lower = xlim[1], upper = xlim[2], resolution = resolution)		
		Dn[i] <- sum(abs(G.n.hat - G.n.tilde))/length(G.n.tilde) #distance measure
	}	
 	x <- which.min(Dn)
	y <- pi0[x]

	if(doplot)
	{
		plot(pi0, Dn, type = 'b') #maybe store in samplesize object
		points(x, y, pch = 19)	
	}
	y
}

##truncate density of effect-sizes for given threshold
truncateEffectsize <- function(obj, threshold)
{
	index.upper <- which.min(abs(obj@theta-threshold))
	lambda <- obj@lambda
	lambda[(length(obj@theta)-index.upper + 1):index.upper] <- 0.0
	lambda/(sum(lambda)*(obj@theta[2]-obj@theta[1]))	
}

##plot density of effect-sizes
plotEffectSize <- function(x, threshold = 0, xlab = "effect size", ylab = "density of effect sizes", main, sub, ...)
{
	if(missing(main)) main <- paste("Experiment:", x@name)

	if(missing(sub)) 
	{
		if(x@method=="Ruppert") # use semi compromise
			sub <- paste("pi0: ", signif(x@pi0[[2]], 3), " (", names(x@pi0[2]),")", sep="")
		else {                #use adjusted pi0 if calculated
			Gm <- ifelse(length(x@pi0) == 1, x@pi0[[1]], x@pi0[[2]])
			sub <- paste("pi0: ", signif(Gm, 3), " (", x@method,")", sep="")
		}
	}

	if(threshold!=0)
		lambda <- truncateEffectsize(x, threshold)
	else
		lambda <- x@lambda

	plot(x@theta, lambda, xlab = xlab, ylab = ylab, main = main, sub = sub, ...)
}

##power estimation
Power <- function(x, threshold = 0, fdr = 0.1, samplesizes = NULL, plot=FALSE, type = 'l', ylim = c(0,1), xlim = c(0,1), 
                  xlab = "p-value", ylab = "average power", main, sub, ...)
{	
	resolution <- 1000 #default
	u <- as.matrix(seq(from = xlim[1], to = xlim[2], length = resolution))	

	if(x@method=="Ruppert")
		Gm <- x@pi0[[2]] # use semi compromise
	else
	  Gm <- ifelse(length(x@pi0) == 1, x@pi0[[1]], x@pi0[[2]]) 	#use adjusted pi0 if calculated

	#Note by Ferreira 
	if(all(fdr >= Gm))
		warning("fdr should be smaller than the proportion of non-differentially expressed genes!")

	r.hat <- fdr*(1 - Gm)/(Gm*(1 - fdr))

	if(is.null(samplesizes))
		sampleSizes <- x@sampleSize 
	else		
		sampleSizes <- samplesizes/2 #effective sample size (1/na + 1/nb)^(-1) assumed na=nb thus na/2

	G.n.hat <- matrix(numeric(), nrow = resolution, ncol = length(sampleSizes))		
	for(i in seq(along = sampleSizes))
		G.n.hat[,i] <- Gnhat(x, sampleSizes[i], threshold=threshold, lower = xlim[1], upper = xlim[2], resolution = resolution)

  if(missing(main)) main <- paste("Experiment:", x@name)
 
	#maybe different way to estimate the intersect
  power <- uStar <- matrix(numeric(length(sampleSizes)*length(fdr)), nrow=length(sampleSizes))
	for(i in 1:length(sampleSizes))
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
	
	#plot the power
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
		legend("bottomright", c(paste("fdr", fdr), paste("sample size", round(2*sampleSizes, 2))), bty='n', lty=lty.legend, pch=pch.legend)

	}	
	colnames(power) <- paste("fdr =", fdr)
	rownames(power) <- paste("sample size =", round(2*sampleSizes, 2))

	power
}


