##' Binning and trimming of the test statistics vector.
##'
##' details follow
##' @title Binning and trimming of the test statistics vector.
##' @param statistics vector of test statistics.
##' @param nbins number of bins.
##' @param trim vector of size two with lower and upper limits for the trimming.
##' @param bin "epdf" or "ecdf" binning using respectively, empirical density of cumulutive distribution function.
##' @param symmetric symmetric trimming TRUE/FALSE.
##' @param plot plot results TRUE/FALSE.
##' @return trimmed and binned vector of test statistics.
##' @author Maarten van Iterson
trimbin <- function(statistics, nbins=100, trim=c(0.01, 0.99), bin=c("epdf", "ecdf"), symmetric=TRUE, plot=TRUE)
  {
    ##trimming  
    q <- quantile(statistics, prob=trim)
    trimmed <- statistics
    
    ##symmetric trimmming
    ##doesn't make sense if quantiles are different
    if(symmetric==TRUE & all.equal(trim[1], 1 - trim[2]) == TRUE)
      {
        q[2] <- max(abs(q[1]), q[2])
        q[1] <- -q[2]
      }
        
    trimmed[statistics <= q[1] | statistics >= q[2]] <- NA    

    ##binning
    breaks <- seq(q[1], q[2], length = nbins + 1)
    x <- breaks[-c(nbins+1)] + 0.5*(q[2] - q[1])/(nbins) #identical to h$mids
    
    if(bin == "epdf")
      {
        if(plot)
          h <- hist(trimmed, breaks = breaks, plot=plot, freq=FALSE)
        else
          h <- hist(trimmed, breaks = breaks, plot=plot)
        y <- h$density
      }
    else if(bin == "ecdf") 
      {
        Fn <- ecdf(trimmed)
        y <- Fn(x)
        if(plot)
          plot(Fn)
      }

    ##define range of the noncentrality parameter in some way?
    list(y=y, x=x)    
  }


#####################################################################################
##Tests:
#####################################################################################

if(FALSE)
{
  statistics <- rnorm(1000)
  tmp <- bintrim(statistics, nbins=100, trim=c(0.01, 0.99), bin="epdf", symmetric=TRUE, plot=TRUE)
  tmp <- bintrim(statistics, nbins=100, trim=c(0.01, 0.99), bin="ecdf", symmetric=TRUE, plot=TRUE)

  statistics <- rchisq(1000, df=3)
  tmp <- bintrim(statistics, nbins=100, trim=c(0.01, 0.90), bin="epdf", symmetric=TRUE, plot=TRUE)
  tmp <- bintrim(statistics, nbins=100, trim=c(0.01, 0.90), bin="ecdf", symmetric=TRUE, plot=TRUE)
  plot(tmp)
}

