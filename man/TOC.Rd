\name{TOC}

\docType{methods}

\alias{TOC}
\alias{TOC,numeric,numeric-method}
\alias{TOC,RasterLayer,RasterLayer-method}

\title{
Construct the table for the TOC curve
}
\description{
Construct the table for the Total Operating Characteristic (TOC) curve for spatial or non-spatial data
}
\usage{
\S4method{TOC}{numeric,numeric}(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, 
progress=FALSE, units=NULL)
\S4method{TOC}{RasterLayer,RasterLayer}(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, 
    progress=FALSE, units=NULL)
}

\arguments{
  \item{index}{
index object of class numeric or RasterLayer
}
  \item{boolean}{
boolean object of class numeric or RasterLayer
}
  \item{mask}{
mask object of class numeric or RasterLayer
}
  \item{nthres}{
an optional integer indicating the number of equal-interval thresholds to be evaluated for the TOC curve. See Details below
}
  \item{thres}{
an optional numeric vector of thresholds to be evaluated for the TOC curve. See Details below
}
  \item{NAval}{
value for nodata (NA values) in the mask object
}
  \item{P}{
count of reference presence observations in the population
}
  \item{Q}{
count of reference absence observations in the population
}
  \item{progress}{
logical; if TRUE, a progress bar is shown
}
  \item{units}{
character string indicating data units
}
}
\details{
thresholds are calculated as the unique values of the index object after masking out NA values (default option), if neither nthres nor thres is provided. The default option can be time-consuming if the amount of unique values in the index object (after masking out NA values) is large (e.g., > 1000). In the latter case, the user may prefer to enter specified thresholds (with the \code{thres} argument), or to indicate the number of equal-interval thresholds to be evaluated for the TOC curve (with the \code{nthres} argument)
}
\value{
an object of class \code{Toc} containing the TOC table, the area under the curve (AUC), maximum AUC and minimum AUC, the prevalence, the population and the data units (for data in the TOC table slot, and the prevalence and population slots) 
}
\seealso{
\code{\link{plot}}
}
\examples{
\dontrun{
index <- raster("p_built01_suitability_1.rst")
boolean <- raster("BuiltGain1985_1999.rst")
mask <- raster("1985NonBuilt01.rst")

# all unique values of the index map after applying the mask are used as thresholds (default 
# option)
tocd <- TOC(index, boolean, mask)
tocd

# thresholds can also be defined by indicating the number of equal-interval thresholds 
tocd <- TOC(index, boolean, mask, nthres=10)
tocd

# A vector of thresholds can be used to define the thresholds
tocd <- TOC(index, boolean, mask, thres=seq(0, 100, by=10))
tocd
}
\dontrun{
index <- raster("Prob_Map2.rst")
boolean <- raster("Change_Map2b.rst")
mask <- raster("MASK3.rst")

# thresholds can also be defined by indicating the number of equal-interval thresholds 
tocd <- TOC(index, boolean, mask, nthres=100)
tocd

# A vector of thresholds can be used to define the thresholds
tocd <- TOC(index, boolean, mask, thres=seq(0, 100, by=10))
tocd
}
}

\keyword{methods}
\keyword{ spatial }
