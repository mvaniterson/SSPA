#######################################################################################
##Class PilotData
#######################################################################################
setClass("PilotData",
	representation(name           = "character",
								 testStatistics = "numeric",                  
                 sampleSize     = "numeric", 	 
								 pValues        = "numeric",
								 dof            = "numeric",
								 nullDist       = "character"),
	prototype(name           = character(1),
            testStatistics = numeric(1),                  
            sampleSize     = numeric(1),       
						pValues        = numeric(1),
  					dof            = numeric(1),
            nullDist       = character(1))
)

##User friendly interface to PilotData
pilotData <- function(name = "Unknown Experiment", 
                      testStatistics = double(1),
                      sampleSizeA = double(1),
                      sampleSizeB = double(1),
											dof = double(1),													
											nullDist = c("normal", "student"))
{
	
	obj <- new("PilotData", name = name)

	if(missing(testStatistics)) 
		stop("Missing test-statistics!")
	if(mode(testStatistics) != mode(numeric(3))) 
		stop("Test statistics not as numeric vector!")
	if(any(is.na(testStatistics)))
		stop("Test statistics contains missing values, not allowed!")
	obj@testStatistics <- testStatistics
	
	#calculate effective sample size
	if(missing(sampleSizeA) && missing(sampleSizeB)) {	
		return(stop("Missing sample sizes!"))
	} else if(missing(sampleSizeA)) {	    
		warning("Assumed equal sample size!")
		obj@sampleSize <- 1/(1/sampleSizeB + 1/sampleSizeB) #sampleSizeB/2
	} else if(missing(sampleSizeB)) {	    
		warning("Assumed equal sample size!")
		obj@sampleSize <- 1/(1/sampleSizeA + 1/sampleSizeA) #sampleSizeA/2
	} else { 
		obj@sampleSize <- 1/(1/sampleSizeA + 1/sampleSizeB)
	}

	obj@nullDist <- match.arg(nullDist)

	if(obj@nullDist == "student" & dof > 0)
	{
		obj@dof <- dof

		#calculate two sided p-values Student
		obj@pValues <- 2*(1 - pt(abs(obj@testStatistics), df=obj@dof, ncp=0))
	} else {	
		#calculate two sided p-values normal
		obj@pValues <- 2*(1 - pnorm(abs(obj@testStatistics), 0, 1))
	}

	obj
}

##Show method for PilotData
setMethod("show", signature("PilotData"), 
function(object) {
	cat("An object of class \"", class(object), "\"\n", sep = "")	
	cat("Experiment name:            ", object@name, "\n", sep = "")
  cat("Number of test-statistics:  ", length(object@testStatistics),"\n", sep = "")
  cat("Effective sample size:      ", round(object@sampleSize, 2),"\n", sep = "")
	cat("Null distribution:      ", object@nullDist,"\n", sep = "")    
	if(object@nullDist=="student")
		cat("Degree of Freedom:      ", object@dof,"\n", sep = "")  
}) 

##Overload Generic hist method for PilotData
setMethod("hist", signature("PilotData"), 
function(x, freq = FALSE, n = 100, xlab = "p-value", main = "Histogram of p-values", ...) { 

	hist(x@pValues, freq = freq, n = n, xlab = xlab,  main = main, ...)
})

##Overload Generic plot method for PilotData (empirical cumulative distribution of p-values)
setMethod("plot", signature(x="PilotData", y="missing"), 
function(x, y, pch=".", col = "black", xlab = "p-value", ylab = "Cumulative Probability", main = "empirical cumulative distribution of p-values",...) {

	plot(sort(x@pValues), (1:length(x@pValues))/length(x@pValues), pch=pch, col = col, xlab = xlab, ylab=ylab, main = main, ...)	
	abline(a=0, b=1, col="black", lty=2)

})

