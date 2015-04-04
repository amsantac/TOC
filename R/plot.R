# Author: Ali Santacruz
# Date :  March 2015
# Version 1.0
# Licence GPL v3


if (!isGeneric("plot")) {
	setGeneric("plot", function(x, ...)
		standardGeneric("plot"))
}	

setMethod("plot", signature='Roc', 
          definition = function(x, labelThres=FALSE, digits=3, modelLeg="Model", ...)  {
            .plotROC(x, labelThres=labelThres, digits=digits, modelLeg=modelLeg, ...)
            return(invisible(NULL))
          }
)


setMethod("plot", signature ='Toc', 
          definition = function(x, labelThres=FALSE, digits=3, modelLeg="Model", ...)  {
            .plotTOC(x, labelThres=labelThres, digits=digits, modelLeg=modelLeg, ...)
            return(invisible(NULL))
          }
)	