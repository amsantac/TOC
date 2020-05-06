\name{ROC}

\docType{methods}

\alias{ROC}
\alias{ROC,numeric,numeric-method}
\alias{ROC,RasterLayer,RasterLayer-method}

\title{
Construct the table for the ROC curve
}
\description{
Construct the table for the Relative Operating Characteristic (ROC) curve for spatial or non-spatial data
}
\usage{
\S4method{ROC}{numeric,numeric}(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, progress=FALSE)

\S4method{ROC}{RasterLayer,RasterLayer}(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, progress=FALSE)
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
an optional integer indicating the number of equal-interval thresholds to be evaluated for the ROC curve. See Details below
}
  \item{thres}{
an optional numeric vector of thresholds to be evaluated for the ROC curve. See Details below
}
  \item{NAval}{
value for nodata (NA values) in the mask object
}
  \item{progress}{
logical; if TRUE, a progress bar is shown
}
}
\details{
thresholds are calculated as the unique values of the index object after masking out NA values (default option), if neither nthres nor thres is provided. The default option can be time-consuming if the amount of unique values in the index object (after masking out NA values) is large (e.g., greater than 1000). In the latter case, the user may prefer to enter specified thresholds (with the \code{thres} argument), or to indicate the number of equal-interval thresholds to be evaluated for the ROC curve (with the \code{nthres} argument)
}
\value{
an object of class \code{Roc} containing the ROC table, the area under the curve (AUC), maximum AUC and minimum AUC
}
\references{
Pontius Jr., R.G., Kangpin, Si. 2014. \emph{The total operating characteristic to measure diagnostic ability for multiple thresholds}. International Journal of Geographical Information Science 28 (3): 570-583. <doi:10.1080/13658816.2013.862623>

Pontius, G., Parmentier, B. 2014. \emph{Recommendations for using the Relative Operating Characteristic (ROC)}. Landscape Ecology 29 (3): 367-382. <doi:10.1007/s10980-013-9984-8>
}
\seealso{
\code{\link{plot}}
}
\examples{
index <- raster(system.file("external/Prob_Map2.rst", package = "TOC"))
boolean <- raster(system.file("external/Change_Map2b.rst", package = "TOC"))
mask <- raster(system.file("external/MASK4.rst", package = "TOC"))

## thresholds can be defined by indicating the number of equal-interval thresholds 
rocd <- ROC(index, boolean, mask, nthres = 100)
rocd

## a vector of thresholds can also be used to define the thresholds
thresholds <- seq(min(unique(index)), max(unique(index)) + 1, 
                  by = ceiling(max(unique(index))/10))
rocd <- ROC(index, boolean, mask, thres = thresholds)
rocd

## all the unique values of the index object can be evaluated as thresholds (default option)
\dontrun{
rocd <- ROC(index, boolean, mask, progress = TRUE)
rocd
}

## generate the ROC curve using non-spatial data (i.e., an object of class numeric)
\dontrun{
index <- getValues(index)
boolean <- getValues(boolean)
mask <- getValues(mask)
rocd <- ROC(index, boolean, mask, nthres = 100)
}
}

\keyword{methods}
\keyword{ spatial }
