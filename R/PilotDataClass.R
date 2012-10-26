#######################################################################################
##Class PilotData
#
#######################################################################################
setClass("PilotData",
	representation(statistics          = "numeric", 
                 effectivesamplesize = "numeric", 	 
								 pvalues             = "numeric",
								 dof                 = "numeric",
								 distribution        = "character"),	
prototype(statistics          = numeric(1),                  
          effectivesamplesize = numeric(1),       
					pvalues             = numeric(1),
  				dof                 = numeric(1),
          distribution        = character(1))
)

#######################################################################################
##Generic Accessor functions for class "PilotData"
#
#######################################################################################
setGeneric("statistics", function(object) { standardGeneric ("statistics") })
setGeneric("effectivesamplesize", function(object) { standardGeneric ("effectivesamplesize") })
setGeneric("pvalues", function(object) { standardGeneric ("pvalues") })
setGeneric("dof", function(object) { standardGeneric ("dof") })
setGeneric("distribution", function(object) { standardGeneric ("distribution") })

setMethod("statistics","PilotData", function(object){ return(object@statistics) })
setMethod("effectivesamplesize","PilotData", function(object){ return(object@effectivesamplesize) })
setMethod("pvalues","PilotData", function(object){ return(object@pvalues) })
setMethod("dof","PilotData", function(object){ return(object@dof) })
setMethod("distribution","PilotData", function(object){ return(object@distribution) })


#######################################################################################
##Validity check method for class "PilotData"
#
#######################################################################################
setValidity("PilotData", function(object){ 
	if(any(is.na(object@statistics)))
		return("Test statistics contains missing values, not allowed!")	
	if(object@distribution == "Student" & object@dof <= 0)
		return("Invalid degrees of freedom!")
	}
)

#######################################################################################
##Show method for class "PilotData"
#
#######################################################################################
setMethod("show", signature("PilotData"), 
function(object) {
	cat("An object of class \"", class(object), "\"\n", sep = "")	
  cat("Number of test-statistics:  ", length(object@statistics),"\n", sep = "")
  cat("Effective sample size:      ", round(object@effectivesamplesize, 2),"\n", sep = "")
	cat("Null distribution:      ", object@distribution,"\n", sep = "")    
	if(object@distribution == "Student")
		cat("Degree of Freedom:      ", object@dof,"\n", sep = "")  
}) 

#######################################################################################
##Overload Generic plot method for class "PilotData" 
#create some diagnostic plots
#######################################################################################
setMethod("plot", signature(x="PilotData"), definition = function(x, y, ...)
{
	object <- x

	dnull <- switch(distribution(object), "Normal" = function(x) dnorm(x), "Student" = function(x) dt(x, df=dof(object)))
	qnull <- switch(distribution(object), "Normal" = function(x) qnorm(x), "Student" = function(x) qt(x, df=dof(object)))

	plot.line <- trellis.par.get("plot.line")
	add.line <- trellis.par.get("add.line")

	hstat <- histogram(statistics(object), col = plot.line$col, border="white", breaks = "Scott", type = "density",
                     ylab="", xlab = "test statistic", main="",
                     panel = function(x, ...) {
                     panel.histogram(x, ...)
									   #panel.mathdensity(dmath = dnull, args = list(), col=add.line$col, lwd=add.line$lwd, lty=add.line$lty)										 
	})

	hpval <- histogram(pvalues(object), breaks = "Scott", col = plot.line$col, border="white", type = "density",
                     ylab="", xlab = "p-value", main = "")

	xypval <- xyplot(seq(along = pvalues(object))/length(pvalues(object)) ~ sort(pvalues(object)), type=c("l", "g"), 
                   ylab="", xlab="p-value (sorted)", main="",
              		 panel = function(x, y, ...){
                   panel.xyplot(x, y, ...)
                   panel.abline(a=0, b=1, col=add.line$col, lwd=add.line$lwd, lty=add.line$lty)
	})	
	
	qqstat <- qqmath(statistics(object), type=c("p", "g"),
                   distribution = qnull,    
                   ylab="test statistic", xlab = ifelse(distribution(object) == "Normal", "qnorm", "qt"), main="",
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

#######################################################################################
##User friendly interface to class "PilotData"
#
#######################################################################################
pilotData <- function(statistics = NULL,
                      effectivesamplesize = NULL,
											dof = double(1),													
											distribution = c("Normal", "Student"))
{
	distribution <- match.arg(distribution)	

	if(missing(statistics))
			stop("Test statistics missing!")	

	if(missing(effectivesamplesize))
			stop("Effective samplesize missing!")	

	if(distribution == "Student" & missing(dof))
			stop("Degrees of freedom missing!")	
	
	#create the object
	object <- new("PilotData", statistics = statistics	, effectivesamplesize = effectivesamplesize, 
                             dof = dof, distribution = distribution)

	if(distribution == "Student") {
		#calculate two sided p-values Student
		object@pvalues <- 2*(1 - pt(abs(statistics), df=dof, ncp=0))
	} else if(distribution == "Normal"){	
		#calculate two sided p-values normal
		object@pvalues <- 2*(1 - pnorm(abs(statistics), 0, 1))
	}		
	object
}


