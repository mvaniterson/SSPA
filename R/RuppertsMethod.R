##Penalized B-Spline regression method of Ruppert et al. 
##Ruppert, D. and Nettleton, D. and Hwang, J.T.G.,
##Exploring the information in p-values for the analysis and planning of multiple-test experiments.,
##Biometrics, 2007, 63, 2, 483-95

##penalized B-splines optimization
constrainedOptimization <- function(lambda, y, X, D)
{
	Dmat <- t(X) %*% X + lambda * t(D) %*% D
	dvec <-  t(X) %*% y

	#define constrains
	Amat <- as.matrix(cbind(rep(1, ncol(X)), diag(rep(1, ncol(X)))))
	bvec <- c(1, rep(0, ncol(X)))

	solve.QP(0.5*Dmat, -1*dvec, Amat, bvec, meq=1, factorized=FALSE)$solution
}

##generalized cross-validation
GCV <- function(lambda, y, X, D, B)
{
	a <- constrainedOptimization(lambda, y, X, D)
	Q <- solve(t(X) %*% X + lambda * t(D) %*% D) #matrix inversion
	DF <- sum(diag(Q %*% (t(X) %*% X )))
	s <- sum((y - X %*% a)^2)
	s /(nrow(X) - DF)^2
}

##Penalized B-Spline regression method of Ruppert et al. 
psreg <- function(object,  resolution = 1000)
{	
	if(!require(quadprog)) warning("In order to use Ruppert's method, the package 'quadprog' must be installed!")
        else require(quadprog)          
	if(!require(splines)) warning("In order to use Ruppert's method, the package 'splines' must be installed!")
        else require(splines)
	
	M <- length(pvalues(object))  	
	w <- 1/(resolution+1)
	x <- seq(w, 1-w, by=w)
	y <- .C("massdist", 
	        x = as.double(pvalues(object)), 
	        xmass = as.double(rep(1/M, M)),
	        nx = as.integer(M), 
	        xlo = as.double(min(x)), 
	        xhi = as.double(max(x)), 
	        y = double(resolution),
  	  	ny = as.integer(resolution),
                PACKAGE = "SSPA")$y/w
	
	if(distribution(object) == "Normal")
	{
		#normal case two-sided
		Zr <- outer(x+w/2, theta(object)*sqrt(effectivesamplesize(object)), function(x, y) 1 - pnorm(qnorm(1-x/2) - y) + pnorm(-qnorm(1-x/2) - y))
		Zl <- outer(x-w/2, theta(object)*sqrt(effectivesamplesize(object)), function(x, y) 1 - pnorm(qnorm(1-x/2) - y) + pnorm(-qnorm(1-x/2) - y))
	} else if(distribution(object) == "Student") {	 
		#student case 
		Zr <- outer(x+w/2, theta(object)*sqrt(effectivesamplesize(object)), 
          function(x, y) 1 - pt(qt(1-x/2, df=dof(object)), df=dof(object), ncp=y) + pt(-qt(1-x/2, df=dof(object)), df=dof(object), ncp=y))
		Zl <- outer(x-w/2, theta(object)*sqrt(effectivesamplesize(object)), 
          function(x, y) 1 - pt(qt(1-x/2, df=dof(object)), df=dof(object), ncp=y) + pt(-qt(1-x/2, df=dof(object)), df=dof(object), ncp=y))
	}

	Z <- (Zr - Zl)/w
	
  #create B-spline basis
	knots <- seq(min(theta(object))+1, max(theta(object))-1, length=object@nKnots+1)	
	B <- splineDesign(knots, x=theta(object), ord=object@bDegree+1, outer.ok = TRUE)
	
	#normalize each B-spline
	B <- B/(colSums(B)*(theta(object)[2] - theta(object)[1]))

	X <- cbind(x, Z%*%B*(theta(object)[2] - theta(object)[1]), deparse.level=0)

	#finite difference second derivative
	D <- diff(diag(ncol(X)), differences=2)

	#generalized cross-validation
	lambda <- optimize(f=GCV, interval=c(0, 1e10), y, X, D)$minimum
	
	gcv <- sapply(10^seq(-10, 12, length=100), GCV, y, X, D, B)	
		

	#optimale solution
	a <- constrainedOptimization(lambda, y, X, D)	

	object@pi0 <- list("semiCompromise" = as.double(X[nrow(X),1:2]%*%a[1:2]),  #semi compromise
	                "semiMin" = as.double(X[nrow(X),]%*%a),                 #semi min
	                "semi" = a[1])		                                      #semi  

	#this is not flexible yet
	b <- a[-1]/(1-as.double(X[nrow(X),]%*%a)) #use semi min

	object@lambda <- as.vector(B %*% b)	

	object
}

