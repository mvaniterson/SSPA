setClass("PilotData", 
         representation(statistics   = "numeric",
                        samplesize   = "numeric",
                        pvalues      = "numeric"),
         prototype(statistics   = numeric(1),
                   samplesize   = numeric(1),
                   pvalues      = numeric(1)),
         contains = ("Distribution"),
         validity = function(object){
           if(any(is.na(object@statistics)))
             return("Test statistics contains missing values, not allowed!")           
           else if(object@samplesize < 0)
             return("Sample size should be large then 0!")
           else
             return(TRUE)}
         )

setGeneric("Statistics", function(object) { standardGeneric ("Statistics") })
setGeneric("Samplesize", function(object) { standardGeneric ("Samplesize") })
setGeneric("Pvalues", function(object) { standardGeneric ("Pvalues") })
setGeneric("distribution", function(object) { standardGeneric ("distribution") })

setMethod("Statistics","PilotData", function(object){ return(object@statistics) })
setMethod("Pvalues","PilotData", function(object){ return(object@pvalues)})
setMethod("Samplesize","PilotData", function(object){ return(object@samplesize)})
setMethod("distribution","PilotData", function(object){ return(object@distribution)})

#######################################################################################
##Show method for class "PilotData"
##
#######################################################################################
setMethod("show", signature("PilotData"), function(object) cat(str(object)))

#######################################################################################
##Overload Generic plot method for class "PilotData"
##create some diagnostic plots
##
#######################################################################################
setMethod("plot", signature(x="PilotData"), definition = function(x, y, ...)
          {
            object <- x

            ##extract functions
            qf0 <- qf0(object)
            
            plot.line <- trellis.par.get("plot.line")
            add.line <- trellis.par.get("add.line")

            hstat <- histogram(Statistics(object), col = plot.line$col, border="white", breaks = "Scott", type = "density",
                               ylab="", xlab = "test statistic", main="",
                               panel = function(x, ...) {
                                 panel.histogram(x, ...)
                                 ##panel.mathdensity(dmath = df0(object), col=add.line$col, lwd=add.line$lwd, lty=add.line$lty) ##doesn't work?
                               })

            hpval <- histogram(Pvalues(object), breaks = "Scott", col = plot.line$col, border="white", type = "density",
                               ylab="", xlab = "p-value", main = "")

            xypval <- xyplot(seq(along = Pvalues(object))/length(Pvalues(object)) ~ sort(Pvalues(object)), type=c("l", "g"),
                             ylab="", xlab="p-value (sorted)", main="",
                             panel = function(x, y, ...){
                               panel.xyplot(x, y, ...)
                               panel.abline(a=0, b=1, col=add.line$col, lwd=add.line$lwd, lty=add.line$lty)
                             })

            qqstat <- qqmath(Statistics(object), type=c("p", "g"),
                             distribution = function(x) qf0(x),
                             ylab="test statistic", xlab = paste("q", distribution(object), sep=""), main="",
                             prepanel = prepanel.qqmathline,
                             panel = function(x, ...) {
                               panel.qqmathline(x, col=add.line$col, lwd=add.line$lwd, lty=add.line$lty, ...)
                               panel.qqmath(x, ...)
                             })

            print(hstat, split=c(1, 1, 2, 2), more=TRUE)
            print(qqstat, split=c(1, 2, 2, 2), more=TRUE)
            print(hpval, split=c(2, 1, 2, 2), more=TRUE)
            print(xypval, split=c(2, 2, 2, 2))

            invisible()
          })

##' User friendly interface to class "PilotData"
##'
##' details follow
##' In the two-group case the effective sample size is defined as the square-root of the inverse of 1/n1 + 1/n2.
##' @title User friendly interface to class "PilotData"
##' @param statistics vector of test statistics
##' @param samplesize total sample size of the pilot-data or effective sample size in two-group case (see Details for more information).
##' @param distribution type of the null/alternative distribution, one of 'norm', 't', 'f' or 'chisq'
##' @param ... additional arguments for the distribution like degrees of freedom
##' @return object of class "PilotData"
##' @author Maarten van Iterson
##' @export
##' @examples
##' pd <- pilotData(statistics=rnorm(100), samplesize=10, distribution="norm")
##' pd
##' plot(pd)
pilotData <- function(statistics = NULL, samplesize = NULL, distribution = c("norm", "t", "f", "chisq"), ...)
{

  distribution <- match.arg(distribution)

  args <- list(...)

  if(missing(statistics))
    stop("Test statistics missing!")

  if(missing(samplesize))
    stop("sample size missing!")

  DistObject <- new("Distribution", distribution=distribution, args=args)
  
  ##create the object
  object <- new("PilotData", statistics = unname(statistics), samplesize = unname(samplesize),
                DistObject)

  ##extract function for calculating the p-values from the Distribution object
  calcPval <- pvalue(DistObject)
  
  object@pvalues <- calcPval(unname(statistics))

  object
}


