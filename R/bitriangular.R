##' Density function for a bi-triangular random variable.
##'
##' For more details see M. Langaas et al. JRSS B 2005.
##' @title Density function for a bi-triangular random variable.
##' @param x vector
##' @param a location of point ... Default a = log2(1.2).
##' @param b location of point ... Default b = log2(4).
##' @param m location of the midpoint of the triangle. Default m = log2(2).
##' @return  Gives the density function.
##' @author Maarten van Iterson
##' @examples
##' curve(dbitri, -4, 4)
##' @export
dbitri <- function(x, a=log2(1.2), b=log2(4), m=log2(2))
{
  y <- numeric(length(x))
  idx <- x >= -b & x < -m
  y[idx] <- (x[idx] + b)/((b-a)*(b-m))

  idx <- x >= -m & x <= -a
  y[idx] <- (x[idx] + a)/((b-a)*(a-m))

  idx <- x >= a & x < m
  y[idx] <-  (x[idx] - a)/((b-a)*(m-a))

  idx <- x >= m & x <= b
  y[idx] <- (x[idx] - b)/((b-a)*(m-b))

  y
}

##' Distribution function for a bi-triangular random variable.
##'
##' For more details see M. Langaas et al. JRSS B 2005.
##' @title Distribution function for a bi-triangular random variable.
##' @param q vector of quantiles.
##' @param a location of point, ... Default a = log2(1.2).
##' @param b location of point, ... Default b = log2(4).
##' @param m location of the midpoint of the triangle. Default m = log2(2).
##' @return Gives the distribution function.
##' @author Maarten van Iterson
##' @examples
##' curve(pbitri, -4, 4)
##' @export
pbitri <- function(q, a=log2(1.2), b=log2(4), m=log2(2))
{
  y <- numeric(length(q))
  idx <- q >= -b & q < -m
  y[idx] <- ((q[idx] + b)^2)/(2*(b-a)*(b-m))

  idx <- q >= -m & q <= -a
  y[idx] <- 0.5 + ((q[idx] + a)^2)/(2*(b-a)*(a-m))

  idx <- q >= -a & q <= a
  y[idx] <- 0.5

  idx <- q >= a & q < m
  y[idx] <- 0.5 + ((q[idx] - a)^2)/(2*(b-a)*(m-a))

  idx <- q >= m & q <= b
  y[idx] <- 1 + ((q[idx] - b)^2)/(2*(b-a)*(m-b))

  idx <- q >= b
  y[idx] <- 1

  y
}

##' Quantile function for a bi-triangular random variable.
##'
##' For more details see M. Langaas et al. JRSS B 2005.
##' @title Quantile function for a bi-triangular random variable.
##' @param p vector of probabilities.
##' @param a location of point, ... Default a = log2(1.2).
##' @param b location of point, ... Default b = log2(4).
##' @param m location of the midpoint of the triangle. Default m = log2(2).
##' @return Gives the quantile function.
##' @author Maarten van Iterson
##' @examples
##' curve(qbitri, 0, 1)
##' @export
qbitri <- function(p, a=log2(1.2), b=log2(4), m=log2(2))
{
  y <- numeric(length(p))

  idx <- p >= 0 & p < pbitri(-m, a=a, b=b, m=m)
  y[idx] <- -b + sqrt(p[idx]*2*(b-a)*(b-m))

  idx <- p >= pbitri(-m, a=a, b=b, m=m) & p < 0.5
  y[idx] <- -a - sqrt((p[idx]- 0.5)*2*(b-a)*(a-m))

  idx <- p == 0.5
  if(sum(idx) > 0)
    y[idx] <- runif(sum(idx), -a, a)

  idx <- p > 0.5 & p <= pbitri(m, a=a, b=b, m=m)
  y[idx] <-  a + sqrt((p[idx] - 0.5)*2*(b-a)*(m-a))

  idx <- p > pbitri(m, a=a, b=b, m=m) & p < 1
  y[idx] <- b - sqrt((p[idx] - 1)*2*(b-a)*(m-b))

  idx <- p >= 1
  y[idx] <- b

  y
}

##' Random generation of bitriangular distributed values.
##'
##' For more details see M. Langaas et al. JRSS B 2005.
##' @title Random generation of bitriangular distributed values.
##' @param n number of observations.
##' @param a location of point, ... Default a = log2(1.2).
##' @param b location of point, ... Default b = log2(4).
##' @param m location of the midpoint of the triangle. Default m = log2(2).
##' @return Generates random deviates.
##' @author Maarten van Iterson
##' @examples
##' hist(rbitri(100), freq=FALSE)
##' curve(dbitri, add=TRUE)
##' @export
rbitri <- function(n, a=log2(1.2), b=log2(4), m=log2(2)) qbitri(runif(n), a=a, b=b, m=m)

##' Simulated microarray data.
##'
##' details follow
##' @title Generate simulated microarray data using the bitriangular distribution.
##' @param mu vector of effect sizes drawn from the bitriangular distribution.
##' @param m number of features (genes, tags, ...).
##' @param pi0 proportion of nondifferentially expressed features.
##' @param J number of samples per group.
##' @param nullX the distribution of nondifferentially expressed features.
##' @param nullY the distribution of nondifferentially expressed features.
##' @param noise standard deviation of the additive noise.
##' @return Matrix of size m x (2J), containing the simulated values.
##' @author Maarten van Iterson
##' @examples
##' ##generate two-group microarray data
##' m <- 5000 ##number of genes
##' J <- 10 ##sample size per group
##' pi0 <- 0.8 ##proportion of non-differentially expressed genes
##' m0 <- as.integer(m*pi0)
##' mu <- rbitri(m - m0, a = log2(1.2), b = log2(4), m = log2(2)) #effect size distribution
##' data <- simdat(mu, m=m, pi0=pi0, J=J, noise=0.01)
##' @export
simdat <- function(mu, m, pi0, J, nullX=function(x)rnorm(x, 0, 1), nullY=function(x)rnorm(x, 0, 1), noise=0.01)
{
  m0 <- as.integer(m*pi0)            #number of true null genes
  m1 <- m - m0                       #number of true alternative genes

  ##Generate Data
  Xm0 <- matrix(nullX(m0*J), ncol=J)
  Ym0 <- matrix(nullY(m0*J), ncol=J)

  Xm1 <- matrix(numeric(m1*J), ncol=J)
  Ym1 <- matrix(numeric(m1*J), ncol=J)
  for(k in 1:m1)
    {
      Xm1[k, ] <- mu[k]/2 + nullX(J)
      Ym1[k, ] <- -1*mu[k]/2 + nullY(J)
    }

  X <- rbind(Xm0, Xm1)
  Y <- rbind(Ym0, Ym1)

  ##error model
  if(!is.null(noise))
    {
      X <- log(exp(rnorm(m, 0, noise)) + exp(X))
      Y <- log(exp(rnorm(m, 0, noise)) + exp(Y))
    }

  cbind(Y=Y, X=X)
}

################################################################################
##some tests
################################################################################
show <- FALSE
if(show)
{
  curve(dbitri, -3, 3, n=1000)
  curve(dbitri(x, b=3, m=1.5, a=log2(1.2)), -4, 4, n=1000)

  integrate(dbitri, -4, 4, b=3, m=1.5, a=log2(1.2))

  curve(pbitri, -3, 3, n=1000)
  curve(pbitri(x, b=3, m=1.5, a=log2(1.2)), -4, 4, n=1000)

  x <- seq(-4, 4, 0.1)
  y <- pbitri(x, b=3, m=2, a=log2(1.2))
  plot(y,x, type='b')
  curve(qbitri(x, b=3, m=2, a=log2(1.2)), add=TRUE, col='red', n=1000)

  curve(qbitri(x, b=3, m=1.5, a=log2(1.2)), 0, 1, n=1000)

  hist(rbitri(1000), n=100)
  hist(rbitri(1000, b=3, m=1.5, a=log2(1.2)), n=100)
}
