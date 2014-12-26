\name{plot.TOC}
\alias{plot.TOC}

\title{
generates a TOC plot
}
\description{
generates a TOC plot
}
\usage{
plot.TOC(toc, labelThres=FALSE, digits)
}

\arguments{
  \item{toc}{
an object of class TOC
}
  \item{labelThres}{
logical, default to FALSE. If TRUE, thresholds are labeled in the TOC plot
}
  \item{digits}{
integer indicating the number of decimal places (round) or significant digits (signif) to be used for labeling the thresholds. Negative values are allowed. See Details in the round function
}
}
\value{
a TOC plot
}
\seealso{
\code{\link{TOC}}
}
\examples{
index <- raster(system.file("external/p_built01_suitability_1.rst", package="TOC"))
boolean <- raster(system.file("external/BuiltGain1985_1999.rst", package="TOC"))
mask <- raster(system.file("external/1985NonBuilt01.rst", package="TOC")) 
tocd <- TOC(index, boolean, mask, NAval=0, uncertainty=TRUE)
tocd
plot.TOC(tocd, labelThres=FALSE)

\dontrun{
index <- raster(system.file("external/Prob_Map2.rst", package="TOC"))
boolean <- raster(system.file("external/Change_Map2b.rst", package="TOC"))
mask <- raster(system.file("external/MASK3.rst", package="TOC"))
tocd <- TOC(index, boolean, mask, nthres=10, NAval=0, uncertainty=TRUE)
tocd
plot.TOC(tocd, labelThres=FALSE)
}
}

\keyword{ spatial }
