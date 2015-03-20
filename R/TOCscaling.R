TOCscaling <- function(tocd, scalingFactor, newUnits){
  #if(!inherits(tocd, TOC) stop("tocd must be an object of class TOC"))
  if(!(is.numeric(scalingFactor))) stop("scalingFactor must be numeric")
  if(!(is.character(newUnits))) stop("newUnits must be a string")
  
  tocd$TOCtable[,2:3] <- tocd$TOCtable[,2:3]/scalingFactor
  tocd$prevalence <- tocd$prevalence/scalingFactor
  tocd$population <- tocd$population/scalingFactor
  tocd$units <- newUnits                              
  return(tocd)
}