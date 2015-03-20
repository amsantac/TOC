\name{plot.ROC}
\alias{plot.ROC}

\title{
Plot the ROC curve
}
\description{
Plot the Relative Operating Characteristic (ROC) curve
}
\usage{
plot.ROC(roc, labelThres=FALSE, digits, modelLeg="Model", ...)
}

\arguments{
  \item{roc}{
an object of class ROC
}
  \item{labelThres}{
logical, default to FALSE. If TRUE, thresholds are labeled in the ROC plot
}
  \item{digits}{
integer indicating the number of decimal places (round) or significant digits (signif) to be used for labeling the thresholds. Negative values are allowed. See Details in the round function
}
  \item{modelLeg}{
a string for labeling the model in the legend 
}
  \item{...}{
additional parameters to be passed to \code{plot}
}
}
\value{
a plot showing the ROC curve
}
\seealso{
\code{\link{ROC}}
}
\examples{
index <- raster(system.file("external/p_built01_suitability_1.rst", package="TOC"))
boolean <- raster(system.file("external/BuiltGain1985_1999.rst", package="TOC"))
mask <- raster(system.file("external/1985NonBuilt01.rst", package="TOC")) 
rocd <- ROC(index, boolean, mask, NAval=0, uncertainty=TRUE)
rocd
plot.ROC(rocd, labelThres=FALSE)

\dontrun{
index <- raster(system.file("external/Prob_Map2.rst", package="TOC"))
boolean <- raster(system.file("external/Change_Map2b.rst", package="TOC"))
mask <- raster(system.file("external/MASK3.rst", package="TOC"))
rocd <- ROC(index, boolean, mask, nthres=100, NAval=0, uncertainty=TRUE)
rocd
plot.ROC(rocd, labelThres=FALSE)
}
}

\keyword{ spatial }
