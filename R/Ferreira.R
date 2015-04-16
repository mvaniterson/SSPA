Gnhat <- function(object, samplesize, threshold = 0, lower = 0, upper = 1.0, resolution)
{
  u <- seq(from = lower + .Machine$double.eps, to = upper - .Machine$double.eps, length = resolution)

  if(missing(samplesize)) samplesize <- Samplesize(object)

  if(threshold != 0)
    lambda <- truncateEffectsize(object, threshold)
  else
    lambda <- Lambda(object)

  theta <- Theta(object)
  ##theta is symmetric thus time can be reduced by a half
  ##theta <- seq(0, max(theta(object)), length=length(theta(object))/2)

  if(distribution(object)=="Normal")
    {
      ##normal case
      gnhat <- outer(qnorm(1-u/2, 0, 1), theta*samplesize, function(x, y) 1 - pnorm(x-y) + pnorm(-x-y))
    } else if(distribution(object)=="Student"){
      ##student case THIS NOT "the shifted-t" but the exact noncentral t

      ## and using the t-dist is slow
      #gnhat <- outer(qt(1-u/2, df=dof(object), ncp=0), theta*sqrt(samplesize),
      #               function(x, y) 1 - pt(x, df=dof(object), ncp=y) + pt(-x, df=dof(object), ncp=y))
    }

  ##theta is symmetric thus time can be reduced by a half
  ##gnhat <- cbind(gnhat[,length(theta):1], gnhat[,1:length(theta)])

  gnhat%*%(lambda*(theta[2]-theta[1]))
}

Hntilde <- function(pvalues, resolution)
{
  n <- length(pvalues)
  x <- seq(from = 0, to = 1, length=resolution)
  m <- round(n^(2/3))
  edf <- sapply(1:m, function(x) length(pvalues[pvalues <= x/m])/n, USE.NAMES = FALSE)
  return(rowSums(matrix(edf*dbinom(1:m, m, rep(x, each = m)), ncol = m, byrow = TRUE)))
}

Gntilde <- function(H.n.tilde, pi0)
{
  F <- seq(from = 0, to = 1, length = length(H.n.tilde)) #uniform distribution
  return((H.n.tilde - F*pi0)/(1 - pi0))
}

Dn <- function(object, pi0, xlim = c(0,1), resolution = 1000, doplot = FALSE)
{
  H.n.tilde <- Hntilde(Pvalues(object), resolution)
  Dn <- numeric(length(pi0))
  for(i in seq(along = pi0))
    {
      object@pi0[["Ferreira"]] <- pi0[i]
      x <- deconvolution(object)
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
