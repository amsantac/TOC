\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,Toc-method}
\alias{plot,Roc-method}

\title{Plot an object of class Toc or Roc}

\description{
Plot a Total Operating Characteristic (TOC) curve or a Relative Operating Characteristic (ROC)
}

\usage{
\S4method{plot}{Toc}(x, labelThres=FALSE, digits=3, modelLeg="Model", ...)
   

\S4method{plot}{Roc}(x, labelThres=FALSE, digits=3, modelLeg="Model", ...)
}

\arguments{
  \item{x}{An object of class Toc or Roc}
  \item{labelThres}{
logical, default to FALSE. If TRUE, thresholds are labeled in the TOC plot
}
  \item{digits}{
integer indicating the number of decimal places (round) or significant digits (signif) to be used for labeling the thresholds and the numeric axes. Negative values are allowed. See Details in the \code{round} function
}
  \item{modelLeg}{
a character string for labeling the model in the legend 
}
\item{...}{
additional parameters to be passed to \code{plot}
}
}

\value{
a plot showing the TOC or the ROC curves
}
\seealso{
\code{\link{TOC}}, \code{\link{ROC}}
}
\examples{
\dontrun{ 
index <- raster("p_built01_suitability_1.rst")
boolean <- raster("BuiltGain1985_1999.rst")
mask <- raster("1985NonBuilt01.rst")
tocd <- TOC(index, boolean, mask)
plot(tocd, main="TOC curve")
rocd <- ROC(index, boolean, mask)
plot(rocd, main="ROC curve")
}
\dontrun{ 
index <- raster("Prob_Map2.rst")
boolean <- raster("Change_Map2b.rst")
mask <- raster("MASK3.rst")
tocd <- TOC(index, boolean, mask, nthres=100)
plot(tocd, main="TOC curve")
rocd <- ROC(index, boolean, mask, nthres=100)
plot(rocd, main="ROC curve")
}
}

\keyword{methods}
\keyword{spatial}

