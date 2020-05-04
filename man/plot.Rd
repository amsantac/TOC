\name{plot}

\docType{methods}

\alias{plot}
\alias{plot,Toc-method}
\alias{plot,Roc-method}

\title{Plot an object of class Toc or Roc}

\description{
Plot a Total Operating Characteristic (TOC) curve or a Relative Operating Characteristic (ROC) curve
}

\usage{
\S4method{plot}{Toc}(x, labelThres=FALSE, modelLeg="Model", digits=3, nticks=5, digitsL=1, 
posL = NULL, offsetL = 0.5, ...)
   

\S4method{plot}{Roc}(x, labelThres=FALSE, modelLeg="Model", digits=3, nticks=5, digitsL=1, 
posL = NULL, offsetL = 0.5, ...)
}

\arguments{
  \item{x}{An object of class Toc or Roc}
  \item{labelThres}{
logical, default to FALSE. If TRUE, thresholds are labeled in the TOC plot
}
  \item{modelLeg}{
a character string for labeling the model in the legend 
}
  \item{digits}{
integer indicating the number of decimal places (round) or significant digits (signif) to be used for labeling the numeric axes. Negative values are allowed. See Details in the \code{round} function
}
  \item{nticks}{
number of tickmarcks to be drawn along the axes
}
  \item{digitsL}{
integer indicating the number of decimal places (round) or significant digits (signif) to be used for labeling the thresholds. Negative values are allowed. See Details in the \code{round} function
}
  \item{posL}{ 
a position specifier for the text labels. Values of 1, 2, 3 and 4, respectively indicate positions below, to the left of, above and to the right of the corresponding coordinates
}
  \item{offsetL}{
when posL is specified, this value gives the offset of the label from the corresponding coordinate in fractions of a character width
}
\item{...}{
additional parameters to be passed to \code{plot}, \code{axis} or \code{text} 
}
}

\value{
a plot showing the TOC or the ROC curve
}
\references{
Pontius Jr., R.G., Kangpin, Si. 2014. \emph{The total operating characteristic to measure diagnostic ability for multiple thresholds}. International Journal of Geographical Information Science 28 (3): 570-583.

Pontius, G., Parmentier, B. 2014. \emph{Recommendations for using the Relative Operating Characteristic (ROC)}. Landscape Ecology 29 (3): 367-382.
}
\seealso{
\code{\link{TOC}}, \code{\link{ROC}}
}
\examples{
index <- raster(system.file("external/Prob_Map2.rst", package = "TOC"))
boolean <- raster(system.file("external/Change_Map2b.rst", package = "TOC"))
mask <- raster(system.file("external/MASK4.rst", package="TOC"))

## create and plot the TOC curve
tocd <- TOC(index, boolean, mask, nthres = 100)
plot(tocd, main = "TOC curve")

## create and plot the ROC curve
rocd <- ROC(index, boolean, mask, nthres = 100)
plot(rocd, main = "ROC curve")

## label the thresholds in the plot
tocd <- TOC(index, boolean, mask, nthres = 10)
plot(tocd, labelThres = TRUE, cex = 0.8, posL = 4)
}

\keyword{methods}
\keyword{spatial}

