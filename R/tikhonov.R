##' Transforms an objective function containing non-diagonal penalty to standard form (diagonal penalty)
##'
##' details follow.
##' @title Transform to standard from
##' @param b the b vector of the system: Ax = b.
##' @param A the A matrix of the system: Ax = b.
##' @param L non-diagonal penalty
##' @return Transformed system ||Ax - b || + lambda ||L|| with L diagonal.
##' @author Maarten van Iterson
standardform <- function(b, A, L)
  {   
    n <- ncol(L);  p <- nrow(L)    
    QR1 <- qr(t(L))
    Kn <- qr.Q(QR1, complete=TRUE) #nxn
    Kp <- qr.Q(QR1, complete=FALSE) #nxp
    Ko <- Kn[, -c(1:p), drop=FALSE]
    Rp <- qr.R(QR1, complete=FALSE) #pxp
    
    piL <- Kp%*%t(solve(Rp))
    
    QR2 <- qr(A%*%Ko)
    Hn <- qr.Q(QR2, complete=TRUE) 
    Ho <- qr.Q(QR2, complete=FALSE)
    Hq <- Hn[, -c(1:(n-p)), drop=FALSE]    
    To <- qr.R(QR2, complete=FALSE) 

    As <- t(Hq)%*%A%*%piL
    bs <- t(Hq)%*%b

    ##back transformation of xs
    ##piL %*% xs + Ko%*%solve(To)%*%t(Ho)%*%(b - A%*%piL%*%xs)
      
    list(As = As, bs = bs, T1 = piL, T2 = Ko%*%solve(To)%*%t(Ho))
  }

##' Tikhonov regularization.
##'
##' details follow.
##' @title Tikhonov regularization
##' @param b the b vector of the system: Ax = b. 
##' @param A the A matrix of the system: Ax = b.
##' @param lambda grid of values for the penalty.
##' @param penalty penalty either 0 = ridge, 1 = first order differences or 2 = second order differences.
##' @return regression coefficients, effective degrees of freedom, intermediate matrix for L-curvature calculation.
##' @author Maarten van Iterson
tikhonov <- function(b, A, lambda, penalty=0)
  {
    ##Ao <- A
    ##bo <- b
    L <- diag(nrow(A))
      
    if(penalty > 0)
      L[1,1] <- 0
    
    ##first or second order difference penalty
    ##if(penalty > 0)
    ##  {
    ##    L <- diff(L, diff=penalty)
    ##    sf <- standardform(b, A, L)
    ##    A <- sf$As
    ##    b <- sf$b
    ##  }
  
    SVD <- svd(A)
    U <- SVD$u; V <- SVD$v; d <- SVD$d 
  
    beta <- z <- matrix(0, nrow=ncol(A), ncol=length(lambda))
    edf <- numeric(length(lambda))
    for(i in 1:length(lambda))
      {
        if(penalty == 0)
          {
            ##svd approach
            beta[,i] <- V %*% diag(d/(d^2 + lambda[i]^2)) %*% t(U) %*% b
            edf[i] <- sum(1/(1 + lambda[i]^2/d^2)) #effective degrees of freedom       
            z[,i] <- V %*% diag(d/(d^2 + lambda[i]^2)) %*% t(U) %*% (A%*%beta[,i] -  b) # for finding the corner of the L-curve
          }
        else
          {
          
            ##inverse explicit probably singular
            tmp <- solve(crossprod(A) + lambda[i]*crossprod(L))
            beta[,i] <-  tmp%*% t(A) %*% b
            edf[i] <- sum(diag(A %*% tmp %*% t(A)))
            z[,i] <-  tmp %*% t(A)%*%(A%*%beta[,i] - b)
            ##augA <- rbind(A, lambda[i]*diag(ncol(A)))
            ##augb <- rbind(b, numeric(ncol(A)))
            ##beta[,i] <- solve(crossprod(augA)) %*% t(augA)%*%augb
          }
      }    
    
    ##back transformation from standard form
    ##if(penalty > 0)
    ##  {
    ##    T1 <- sf$T1
    ##    T2 <- sf$T2
    ##    beta <- apply(beta, 2, function(x) T1%*%x + T2%*%(bo - Ao%*%T1%*%x))
    ##  }    
    
    list(beta = beta, edf=edf, z=z)
  }

##' Find optimal regularization parameter
##'
##' details follow.
##' @title Find optimal regularization parameter
##' @param b the b vector of the system: Ax = b. 
##' @param A the A matrix of the system: Ax = b.
##' @param beta regression coefficients.
##' @param edf effective degrees of freedom.
##' @param lambda grid of penalty values.
##' @param z intermediate matrix for L-curvature calculation.
##' @param method Either the L-curve, GCV or AIC. 
##' @param plot Plot TRUE/FALSE.
##' @param log Plot on log-scale TRUE/FALSE.
##' @param verbose Verbose TRUE/FALSE
##' @return generates optionally figure and returns the index for the optimal penalty value.
##' @author Maarten van Iterson
regularization <- function(b, A, beta, edf, lambda, z, method=c("lcurve", "gcv", "aic"), plot=TRUE, log=TRUE, verbose=FALSE)
  {   
    ##one liners for regularization    
    normR <- colSums((A%*%beta - b)^2)  #norm of the residuals    
    normS <- colSums(beta^2)            #norm of the solution
    gcv <- normR/(1 - edf/length(b))^2 #generalized cross-validation error
    aic <- log(normR) + 2*edf/length(b) #Akaike's information criterion
    
    ##plot numbers with colors
    col <- (1:length(lambda))%/%10 + 1
    pch <- as.character((1:length(lambda))%%10)
    
    if(method == "gcv")
      {
        idx <- cornerScurve(A, beta, b, lambda, z, verbose=verbose)
        print(paste("Minimum at:", which.min(gcv)))
        if(plot)
          {
            if(log)
              {
                plot(lambda, gcv, log="xy", pch=pch, col=col, ylab=expression(log[10](GCV)), xlab=expression(log[10](lambda)), main="GCV")
                points(lambda[idx], gcv[idx], pch=4, cex=2)                
              }
            else
              {                
                plot(lambda, gcv, pch=pch, col=col, ylab=expression(GCV), xlab=expression(lambda), main="GCV")
                points(lambda[which.min(gcv)], gcv[which.min(gcv)], pch=4, cex=2)
              }
          }
        if(log)
          return(idx)
        else
          return(which.min(gcv))
      }
    if(method == "aic")
      {
        idx <- cornerScurve(A, beta, b, lambda, z, verbose=verbose)
        if(plot)
          {
            plot(lambda, aic, log="x", main="AIC", ylab="AIC", xlab=expression(log[10](lambda)), pch=pch, col=col)
            points(lambda[idx], aic[idx], pch=4, cex=2)
          }      
        return(idx)
      }
    if(method == "lcurve")
      {
        idx <- cornerLcurve(A, beta, b, lambda, z, verbose=verbose)
        if(plot)
          {
            plot(0.5*log10(normR), 0.5*log10(normS), pch=pch, col=col, ylab=expression(log[10](L[2]-norm(beta))),
                 xlab=expression(log[10](L[2]-norm(residuals))), main="L-curve")
            points(0.5*log10(normR)[idx], 0.5*log10(normS)[idx], pch=4, cex=2)
          }      
        return(idx)
      }
    if(method=="new")
      {
        if(plot)
          {
            par(mfcol=c(1,2))  
            plot(lambda, sqrt(edf^2 + normR^2), log="x", pch=pch, col=col, ylab="", xlab="", main="new")
            plot(edf, normR, log="x", pch=pch, col=col, ylab="", xlab="", main="new")
            par(mfcol=c(1,1))
          }
      }
  }

##' Find corner L-curve
##'
##' details follow.
##' @title Find corner L-curve
##' @param A the A matrix of the system: Ax = b.
##' @param beta regression coefficients.
##' @param b the b vector of the system: Ax = b. 
##' @param lambda grid of penalty values.
##' @param z intermediate matrix for L-curvature calculation.
##' @param verbose Verbose TRUE/FALSE
##' @return index for the corner of the L-curve.
##' @author Maarten van Iterson
cornerLcurve <- function(A, beta, b, lambda, z, verbose=FALSE)
  {
    xi <- colSums(beta^2)               #norm of the solution
    rho <- colSums((A%*%beta - b)^2)    #norm of the residuals

    ##dxi 
    dxi <- 4*diag(crossprod(beta, z))/lambda

    ##or dxi
    ##sigma <- svd(A)$d
    ##U <- svd(A)$u
    ##dxi <- numeric(length(lambda))    
    ##Q <- t(b)%*%U
    ##for(i in 1:length(lambda))
    ##  {
    ##    dtmp <- -4*lambda[i]*sigma^2/(sigma^2+lambda[i]^2)^3       
    ##    dxi[i] <- Q%*%diag(dtmp)%*%t(Q)     
    ##}      
    
    ##formula of Hansen Discrete Inverse Problems p 92 contains an error should be lambda^4 in the denominator?
    k <- (2*xi*rho/dxi)*(lambda^2*dxi*rho + 2*lambda*xi*rho + lambda^4*xi*dxi)/(lambda^4*xi^2 + rho^2)^(3/2)
    
    ##original
    ##k <- (2*xi*rho/dxi)*(lambda^2*dxi*rho + 2*lambda*xi*rho + lambda^4*xi*dxi)/(lambda^2*xi^2 + rho^2)^(3/2)

    ##plot numbers with colors
    col <- (1:length(lambda))%/%10 + 1
    pch <- as.character((1:length(lambda))%%10)
   
    ##check curvature
    if(verbose)
      {
        plot(lambda, k, xlab=expression(lambda), ylab="curvature", log="x", type="b", pch=pch, col=col)      
      }
    
    which.min(k)
  }

##' Find corner S-curve
##'
##' details follow.
##' @title Find corner S-curve
##' @param A the A matrix of the system: Ax = b.
##' @param beta regression coefficients.
##' @param b the b vector of the system: Ax = b. 
##' @param lambda grid of penalty values.
##' @param z intermediate matrix for L-curvature calculation.
##' @param verbose Verbose TRUE/FALSE
##' @return index for the corner of the S-curve.
##' @author Maarten van Iterson
cornerScurve <- function(A, beta, b, lambda, z, verbose=FALSE)
  {    
    n <- length(b)    
    sigma <- svd(A)$d
    U <- svd(A)$u
                  
    xi <- colSums(beta^2)               #norm of the solution
    rho <- colSums((A%*%beta - b)^2)    #norm of the residuals
    
    dxi <- 4*diag(crossprod(beta, z))/lambda     
    
    B <- numeric(length(lambda))    
    Q <- t(b)%*%U
    for(i in 1:length(lambda))
      {        
        tmp <- sigma^2/(sigma^2 + lambda[i]^2)^4
        B[i] <- Q%*%diag(tmp)%*%t(Q)
    }    

    ##ddxi <- dxi/lambda + 24 * lambda^2 * ddxi
    
    drho <- -lambda^2*dxi
    ##ddrho <- -2*lambda*dxi - lambda^2*ddxi   

    tau <- sapply(lambda, function(x)sum(sigma^2/(sigma^2+x^2)))    
    dtau <- -2*lambda*sapply(lambda, function(x) sum(sigma^2/(sigma^2+x^2)^2))       
    A <- sapply(lambda, function(x) sum(sigma^2/(sigma^2+x^2)^3))
    
    ##ddtau <- 2*sapply(lambda, function(x) sum((sigma^2*(3*x^2 - sigma^2))/(sigma^2+x^2)^3))
    
    numerator <- -4*dxi/rho - 24*lambda^3*B/rho - (lambda^3)*(dxi/rho)^2 + 4*dtau/(n*lambda^2) + 16 *lambda*A/n
    
    denominator <- (1/lambda^2 + (2*dtau/n -lambda^2*dxi/rho)^2)^(3/2)
    
    k <- numerator/denominator

    ##plot numbers with colors
    col <- (1:length(lambda))%/%10 + 1
    pch <- as.character((1:length(lambda))%%10)  
    
    ##check the curvature
    if(verbose)
      plot(lambda, k, xlab=expression(lambda), ylab="curvature", log="x", type="b", pch=pch, col=col)
        
    which.max(k)
  }

##' Generates Picard-plot
##'
##' details follow.
##' @title Picard-plot
##' @param A the A matrix of the system: Ax = b.
##' @param b the b vector of the system: Ax = b. 
##' @param xlim xlim of Picard-plot.
##' @param ylim ylim of Picard-plot.
##' @param main main
##' @param legend legend
##' @return generates Picard-plot.
##' @author Maarten van Iterson
picardplot <- function(A, b, xlim=c(1, length(b)), ylim=NULL, main="Picard-Plot", legend=TRUE)
  {
    SVD <- svd(A)
    U <- SVD$u
    sigma <- SVD$d     
    
    if(is.null(ylim))
      ylim <- c(min(sigma), 1/min(sigma))

    plot(sigma, log="y", pch=1, type="p", ylim=ylim, xlim=xlim, ylab="", xlab="", main=main)
    points(abs(crossprod(U, b)), pch=3, col=3)
    points(abs(crossprod(U, b))/sigma, pch=6, col=6)
    abline(h=1e-16, lty=2, lwd=2, col=1)
    abline(h=1e0, lty=1, lwd=2, col=1)
    
    if(legend)
      legend("center", expression(sigma[i], u[i]^t*b, frac(u[i]^t*b, sigma[i])), pch=c(1,3,6), col=c(1,3,6), bty="n")
  }


Tikhonov <- function(object)
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
    method <- object@control$modelselection
    log <- object@control$log
    penalty <- object@control$penalty
    lambda <- object@control$lambda

      if(scale == "cdfpval")
        statistics <- Pvalues(object)
      else
        statistics <- Statistics(object)

    resolution <- length(theta)

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
                function(x, y) 1 - pf1(q=qf0(p=1 - y/2), y=N*x) + pf1(q=-qf0(p=1 - y/2), y=N*x) ##two-sided t norm
                else
                function(x, y) 1 - pf1(q=qf0(p=1 - y), y=N*x) ##one-side f chisq
                )
     
    A <- midpoint(K, theta[1], theta[resolution], resolution-1, x)
    A <- cbind(K(0, x), t(A))  

    object@info$A <- A
    object@info$b <- b
    
    tkh <- tikhonov(b, A, lambda, penalty = 0)
    
    object@info$beta <- tkh$beta
    object@info$edf <- tkh$edf
    object@info$z <- tkh$z

    idx <- regularization(b = b, A = A, beta = tkh$beta, edf = tkh$edf, lambda = lambda, z = tkh$z, method=method, plot=verbose, log=log, verbose=verbose)

    object@info$idx <- idx   
    
    beta <- tkh$beta[, idx]

    object@pi0 <- beta[1]
    if(beta[1] > 1)
      warning("Estimated pi0 > 1!")

    C <- sum(beta[-1]/(1-beta[1]))*(theta[2]-theta[1])

    object@lambda <- beta[-1]/(C*(1-beta[1]))

    object@theta <- theta
    object
  }
