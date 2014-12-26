\name{TOC}
\alias{TOC}

\title{
calculates TOC table
}
\description{
calculates TOC
}
\usage{
TOC(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, ranking=FALSE, P=NA, Q=NA, uncertainty=FALSE)
}

\arguments{
  \item{index}{
index Raster map
}
  \item{boolean}{
boolean Raster map
}
  \item{mask}{
mask Raster map
}
  \item{nthres}{
an optional integer indicating the number of thresholds to be evaluated for the TOC curve. See Details below
}
  \item{thres}{
an optional numeric vector of thresholds to be evaluated for the TOC curve. See Details below
}
  \item{NAval}{
value for nodata
}
  \item{ranking}{
logical; default to FALSE. If TRUE, cell values are ranked to solve ties
}
  \item{P}{
count of reference presence observations in the population
}
  \item{Q}{
count of reference absence observations in the population
}
  \item{uncertainty}{
logical; default to FALSE. If TRUE, uncertainty in AUC calculation is computed and maximum AUC and minimum AUC, given the uncertainty, are provided in the output. See Details below
}
}
\details{
thresholds are calculated as the unique values of the index map after masking out NA values (default option), if neither nthres nor thres is provided. The default option can be time-consuming if the amount of unique values in the index map (after masking out NA values) is large (e.g., > 1000) 
}
\value{
a list of class TOC containing the TOC table, the area under the curve (AUC), the map coordinate units, the prevalence (as pixel counts), the population (as pixel counts). Maximum AUC and minimum AUC are also provided if uncertainty is set to TRUE
}
\seealso{
\code{\link{plot.TOC}}
}
\examples{
index <- raster(system.file("external/p_built01_suitability_1.rst", package="TOC"))
boolean <- raster(system.file("external/BuiltGain1985_1999.rst", package="TOC"))
mask <- raster(system.file("external/1985NonBuilt01.rst", package="TOC")) 

# all unique values of the index map after applying the mask are used as thresholds (default option)
tocd <- TOC(index, boolean, mask, NAval=0, uncertainty=TRUE)
tocd

\dontrun{
# thresholds can also be defined by indicating the number of equal-interval thresholds 
tocd <- TOC(index, boolean, mask, nthres=10, NAval=0, uncertainty=TRUE)
tocd

# A vector of thresholds can be used to define the thresholds
tocd <- TOC(index, boolean, mask, thres=seq(0, 100, by=10), NAval=0, uncertainty=TRUE)
tocd
}

\dontrun{
index <- raster(system.file("external/Prob_Map2.rst", package="TOC"))
boolean <- raster(system.file("external/Change_Map2b.rst", package="TOC"))
mask <- raster(system.file("external/MASK3.rst", package="TOC"))
tocd <- TOC(index, boolean, mask, nthres=10, NAval=0, uncertainty=TRUE)
}
}

\keyword{ spatial }
