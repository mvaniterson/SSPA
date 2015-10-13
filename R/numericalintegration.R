##' Implementation of the midpoint rule for the numerical integration of uni- and bivariate functions.
##'
##' details follow
##' @title Midpoint rule for numerical integration.
##' @param f bivariate function.
##' @param a lower bound of the grid.
##' @param b upper bound of the grid.
##' @param n grid size.
##' @param ... trick for evaluating the second parameter in case a bivariate functions was supplied.
##' @return vector or matrix of function evaluations use sum to obtain the integrand.
##' @author Maarten van Iterson
midpoint <- function(f, a, b, n, ...)
  {    
    x <- seq(a, b, length=n)
    ##trick for handling two variable functions
    input <- list(...)
    if(length(input) == 0)
      z <- f(x)
    else
      {
        y <- unlist(input)  
        z <- outer(x, y, f)
      }
    z*(x[2]-x[1])    
  }

##' Implementation of the trapezoidal rule for the numerical integration of uni- and bivariate functions.
##'
##' details follow
##' @title Trapezoidal rule for numrical integration.
##' @param f bivariate function.
##' @param a lower bound of the grid.
##' @param b upper bound of the grid.
##' @param n grid size.
##' @param ... trick for evaluating the second parameter in case a bivariate functions was supplied.
##' @return vector or matrix of function evaluations use sum to obtain the integrand.
##' @author Maarten van Iterson
trapezoidal <- function(f, a, b, n, ...)
  {
    x <- seq(a, b, length=n)
    w <- c(1, rep(2, n-2), 1)
    ##trick for handling two variable functions
    input <- list(...)
    if(length(input) == 0)
      z <- f(x)
    else
      {
        y <- unlist(input)  
        z <- outer(x, y, f)
      }
    z*w*(x[2]-x[1])/2
  }

##' Implementation of Simpson's rule for the numerical integration of uni- and bivariate functions.
##'
##' details follow
##' @title Simpson's rule for numrical integration.
##' @param f bivariate function.
##' @param a lower bound of the grid.
##' @param b upper bound of the grid.
##' @param n grid size.
##' @param ... trick for evaluating the second parameter in case a bivariate functions was supplied. 
##' @return  vector or matrix of function evaluations use sum to obtain the integrand.
##' @author Maarten van Iterson
simpson <- function(f, a, b, n = 5, ...)
  {  
    #n should be odd and >= 5
    if(n < 5)
      n <- 5
    if(n%%2 == 0) #n is even add 1
       n <- n + 1
    
    x <- seq(a, b, length=n)
    w <- c(1, rep(c(4, 2), (n-3)/2), 4, 1)

    ##trick for handling two variable functions
    input <- list(...)
    if(length(input) == 0)
      z <- f(x)
    else
      {
        y <- unlist(input)  
        z <- outer(x, y, f)
      }

    z*w*(x[2]-x[1])/3
  }

