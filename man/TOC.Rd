\name{TOC}

\docType{methods}

\alias{TOC}
\alias{TOC,numeric,numeric-method}
\alias{TOC,RasterLayer,RasterLayer-method}

\title{
Construct the table for the TOC curve
}
\description{
Construct the table for the Total Operating Characteristic (TOC) curve for spatial or non-spatial data. The TOC method is a modification of the ROC method which measures the ability of an index variable to diagnose either presence or absence of a characteristic. The diagnosis depends on whether the value of an index variable is above a threshold. Each threshold generates a two-by-two contingency table, which contains four entries: hits (H), misses (M), false alarms (FA), and correct rejections (CR). While ROC shows for each threshold only two ratios, H/(H + M) and FA/(FA + CR), TOC reveals the size of every entry in the contingency table for each threshold (Pontius Jr., R.G., Si, K. 2014. <doi:10.1080/13658816.2013.862623>).
}
\usage{
\S4method{TOC}{numeric,numeric}(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, 
progress=FALSE, units=character(0))
\S4method{TOC}{RasterLayer,RasterLayer}(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, 
    progress=FALSE)
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
thresholds are calculated as the unique values of the index object after masking out NA values (default option), if neither nthres nor thres is provided. The default option can be time-consuming if the amount of unique values in the index object (after masking out NA values) is large (e.g., greater than 1000). In the latter case, the user may prefer to enter specified thresholds (with the \code{thres} argument), or to indicate the number of equal-interval thresholds to be evaluated for the TOC curve (with the \code{nthres} argument)
}
\value{
an object of class \code{Toc} containing the TOC table, the area under the curve (AUC), maximum AUC and minimum AUC, the prevalence, the population and the data units (for data in the TOC table slot, and the prevalence and population slots) 
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
tocd <- TOC(index, boolean, mask, nthres = 100)
tocd

## a vector of thresholds can also be used to define the thresholds
thresholds <- seq(min(unique(index)), max(unique(index)) + 1, 
                  by = ceiling(max(unique(index))/10))
tocd <- TOC(index, boolean, mask, thres = thresholds)
tocd

## all the unique values of the index object can be evaluated as thresholds (default option)
\dontrun{
tocd <- TOC(index, boolean, mask, progress = TRUE)
tocd
}

## generate the TOC curve using non-spatial data (i.e., an object of class numeric)
\dontrun{
index <- getValues(index)
boolean <- getValues(boolean)
mask <- getValues(mask)
tocd <- TOC(index, boolean, mask, nthres = 100)
}
}

\keyword{methods}
\keyword{ spatial }
