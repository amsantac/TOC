\name{scaling}

\docType{methods}

\alias{scaling}
\alias{scaling,Toc-method}

\title{scale the output TOC values and change units}

\description{
scale the 'Hits' and 'Hits+FalseAlarms' values in the TOC output table, as well as the prevalence and population, using a scaling factor. Labels for the modified units in the TOC object are changed to \code{newUnits}
}

\usage{
\S4method{scaling}{Toc}(x, scalingFactor, newUnits)
   
}

\arguments{
  \item{x}{
an object of class Toc
}
  \item{scalingFactor}{
numeric value to scale 'Hits' and 'Hits+FalseAlarms' values in the TOC output table, as well as the prevalence and population
}
  \item{newUnits}{
charater string for the new data units in the TOC object
}
}
\value{
an object of class TOC
}
\seealso{
\code{\link{TOC}}, \code{\link{ROC}}
}
\examples{
index <- raster(system.file("external/Prob_Map2.rst", package = "TOC"))
boolean <- raster(system.file("external/Change_Map2b.rst", package = "TOC"))
mask <- raster(system.file("external/MASK4.rst", package = "TOC"))
tocd <- TOC(index, boolean, mask, nthres = 100)
plot(tocd)

## scale units from square m to square km
tocd_sqkm <- scaling(tocd, scalingFactor = 1000000, newUnits = "square km")
plot(tocd_sqkm) 
}

\keyword{methods}
\keyword{spatial}

