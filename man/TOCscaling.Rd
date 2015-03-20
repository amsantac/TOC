\name{TOCscaling}
\alias{TOCscaling}

\title{
scale the TOC output values and change units 
}
\description{
scale the 'Hits' and 'Hits+FalseAlarms' values in the TOC output table, as well as the prevalence and population, using a scaling factor. Labels for the modified units in the TOC object are changed to \code{newUnits}
}
\usage{
TOCscaling(tocd, scalingFactor, newUnits)
}

\arguments{
  \item{tocd}{
an object of class TOC
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
\code{\link{TOC}}
}
\examples{
index <- raster(system.file("external/p_built01_suitability_1.rst", package="TOC"))
boolean <- raster(system.file("external/BuiltGain1985_1999.rst", package="TOC"))
mask <- raster(system.file("external/1985NonBuilt01.rst", package="TOC")) 
tocd <- TOC(index, boolean, mask, NAval=0, uncertainty=TRUE)
plot.TOC(tocd, labelThres=FALSE)

# scale units from square m to square km
tocd_sqkm <- TOCscaling(tocd, scalingFactor=1000000, newUnits="square km")
plot.TOC(tocd_sqkm) 

\dontrun{
index <- raster(system.file("external/Prob_Map2.rst", package="TOC"))
boolean <- raster(system.file("external/Change_Map2b.rst", package="TOC"))
mask <- raster(system.file("external/MASK3.rst", package="TOC"))

# scale units from square m to thousand square km
tocd_sqkm <- TOCscaling(tocd, scalingFactor=1000000, newUnits="thousand square km")
plot.TOC(tocd_sqkm)
}
}

\keyword{ spatial }
