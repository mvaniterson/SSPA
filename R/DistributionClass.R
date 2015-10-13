#######################################################################################
##Class Distribution
##
#######################################################################################
setClass("Distribution",
         representation(distribution = "character",
                        args         = "list"),
         prototype(distribution = character(1),
                   args         = list())
         )

#######################################################################################
##Validity check method for class "Distribution"
##
#######################################################################################
setValidity("Distribution", function(object){

  distribution <- object@distribution
  args <- object@args

  if(distribution == "t" | distribution == "chisq")
    {
      if(length(names(args)) != 1)
        return("Expected one argument: `df'!")
      else
        {
          if(names(args) != "df")
            return("Expected argument: `df'!")
        }
    }
  else if(distribution == "f")
    {
      if(length(names(args)) != 2)
        return("Expected two arguments: `df1' and `df2'!")
      else
        {
          if(!all(names(args) %in% c("df1", "df2")))
            return("Expected arguments: `df1' and `df2'!")
        }
    }

}
            )

#######################################################################################
##Generic functions for class "Distribution"
##
#######################################################################################
setGeneric("df0", function(object) { standardGeneric ("df0") })
setGeneric("pf0", function(object) { standardGeneric ("pf0") })
setGeneric("qf0", function(object) { standardGeneric ("qf0") })
setGeneric("cf0", function(object) { standardGeneric ("cf0") })

setGeneric("df1", function(object) { standardGeneric ("df1") })
setGeneric("pf1", function(object) { standardGeneric ("pf1") })

setGeneric("pvalue", function(object) { standardGeneric ("pvalue") })

setMethod("df0","Distribution", function(object){
  distribution <- object@distribution
  args <- object@args
  df0 <- switch(object@distribution,
                norm = function(x) dnorm(x, mean=0, sd=1),
                t = function(x, dof = args$df) dt(x, df=dof, ncp=0),
                f = function(x, dof1 = args$df1, dof2 = args$df2) df(x, df1=dof1, df2=dof2, ncp=0),
                chisq = function(x, dof = args$df) dchisq(x, df=dof, ncp=0)
                )
  return(df0)
})

setMethod("pf0","Distribution", function(object){
  distribution <- object@distribution
  args <- object@args
  pf0 <- switch(object@distribution,
                norm = function(q) pnorm(q, mean=0, sd=1),
                t = function(q, dof = args$df) pt(q, df=dof, ncp=0),
                f = function(q, dof1 = args$df1, dof2 = args$df2) pf(q, df1=dof1, df2=dof2, ncp=0),
                chisq = function(q, dof = args$df) pchisq(q, df=dof, ncp=0)
                )
  return(pf0)
})

setMethod("qf0","Distribution", function(object){
  distribution <- object@distribution
  args <- object@args
  qf0 <- switch(object@distribution,
                norm = function(p) qnorm(p, mean=0, sd=1),
                t = function(p, dof = args$df) qt(p, df=dof, ncp=0),
                f = function(p, dof1 = args$df1, dof2 = args$df2) qf(p, df1=dof1, df2=dof2, ncp=0),
                chisq = function(p, dof = args$df) qchisq(p, df=dof, ncp=0)
                )
  return(qf0)
})

setMethod("cf0","Distribution", function(object){
  distribution <- object@distribution
  args <- object@args
  cf0 <- switch(object@distribution,
                norm = function(t) cnorm(t, mean=0, sd=1),
                t = function(t, dof = args$df) ct(t, df=dof),
                f = function(t, dof1 = args$df1, dof2 = args$df2) stop("Characteristic function for the F not implemented!"),
                chisq = function(t, dof = args$df) stop("Characteristic function for the chisq not implemented!"),
                )
  return(cf0)
})

##suppressWarnings to get rid of the annoying "full precision may not have been achieved in 'pnt{final}"
##TODO maybe come up with nicer solution

setMethod("df1","Distribution", function(object){
  distribution <- object@distribution
  args <- object@args
  df1 <- switch(object@distribution,
                norm = function(x, y) dnorm(x, mean=y, sd=1),
                t = function(x, y, dof = args$df) suppressWarnings(dt(x, df=dof, ncp=y)),
                f = function(x, y, dof1 = args$df1, dof2 = args$df2) df(x, df1=dof1, df2=dof2, ncp=y),
                chisq = function(x, y, dof = args$df) dchisq(x, df=dof, ncp=y)
                )
  return(df1)
})

setMethod("pf1","Distribution", function(object){
  distribution <- object@distribution
  args <- object@args
  pf1 <- switch(object@distribution,
                norm = function(q, y) pnorm(q, mean=y, sd=1),
                t = function(q, y, dof = args$df) suppressWarnings(pt(q, df=dof, ncp=y)),
                f = function(q, y, dof1 = args$df1, dof2 = args$df2) pf(q, df1=dof1, df2=dof2, ncp=y),
                chisq = function(q, y, dof = args$df) pchisq(q, df=dof, ncp=y)
                )
  return(pf1)
})

setMethod("pvalue", "Distribution", function(object){
  distribution <- object@distribution
  if(distribution == "norm" | distribution == "t")
    return(function(x) 2*(1 - pf0(object)(abs(x)))) ##two-sided
  else
    return(function(x) 1 - pf0(object)(x)) ##one-sided
})

#######################################################################################
##Show method for class "Distribution"
##
#######################################################################################
setMethod("show", signature("Distribution"),
          function(object) {
            cat("An object of class \"", class(object), "\"\n", sep = "")
            cat("Distribution:  ", object@distribution,"\n", sep = "")
            if(length(object@args) > 0)
              {
                cat("additional arguments:      ")
                for(i in 1:length(object@args))
                  cat(names(object@args)[i], "=", object@args[[1]], "\n")                                              }
          })

