# Author: Ali Santacruz
# Date :  March 2015
# Version 1.0
# Licence GPL v3


# set x@table
if (!isGeneric("setTable")) {
  setGeneric("setTable", function(x, ...)
    standardGeneric("setTable"))
}	

setMethod("setTable", signature = 'Roc', 
          def = function(x, value)  {
            x@table <- value
            validObject(x)
            return(x)
          }
)

setMethod("setTable", signature = 'Toc', 
          def = function(x, value)  {
            x@table <- value
            validObject(x)
            return(x)
          }
)

# set x@AUC
if (!isGeneric("setAUC")) {
  setGeneric("setAUC", function(x, ...)
    standardGeneric("setAUC"))
}  

setMethod("setAUC", signature = 'Roc', 
          def = function(x, value)  {
            x@AUC <- value
            validObject(x)
            return(x)
          }
)

setMethod("setAUC", signature = 'Toc', 
          def = function(x, value)  {
            x@AUC <- value
            validObject(x)
            return(x)
          }
)

# set x@maxAUC
if (!isGeneric("setMaxAUC")) {
  setGeneric("setMaxAUC", function(x, ...)
    standardGeneric("setMaxAUC"))
}  

setMethod("setMaxAUC", signature = 'Roc', 
          def = function(x, value)  {
            x@maxAUC <- value
            validObject(x)
            return(x)
          }
)

setMethod("setMaxAUC", signature = 'Toc', 
          def = function(x, value)  {
            x@maxAUC <- value
            validObject(x)
            return(x)
          }
)

# set x@minAUC
if (!isGeneric("setMinAUC")) {
  setGeneric("setMinAUC", function(x, ...)
    standardGeneric("setMinAUC"))
}  

setMethod("setMinAUC", signature = 'Roc', 
          def = function(x, value)  {
            x@minAUC <- value
            validObject(x)
            return(x)
          }
)

setMethod("setMinAUC", signature = 'Toc', 
          def = function(x, value)  {
            x@minAUC <- value
            validObject(x)
            return(x)
          }
)

# set x@prevalence
if (!isGeneric("setPrevalence")) {
  setGeneric("setPrevalence", function(x, ...)
    standardGeneric("setPrevalence"))
}  

setMethod("setPrevalence", signature = 'Toc', 
          def = function(x, value)  {
            x@prevalence <- value
            validObject(x)
            return(x)
          }
)

# set x@population
if (!isGeneric("setPopulation")) {
  setGeneric("setPopulation", function(x, ...)
    standardGeneric("setPopulation"))
}  

setMethod("setPopulation", signature = 'Toc', 
          def = function(x, value)  {
            x@population <- value
            validObject(x)
            return(x)
          }
)

# set x@units
if (!isGeneric("setUnits")) {
  setGeneric("setUnits", function(x, ...)
    standardGeneric("setUnits"))
}  

setMethod("setUnits", signature = 'Toc', 
          def = function(x, value)  {
            x@units <- value
            validObject(x)
            return(x)
          }
)
