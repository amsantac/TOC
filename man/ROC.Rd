\name{ROC}
\alias{ROC}

\title{
Construct the table for the ROC curve
}
\description{
Construct the table for the Relative Operating Characteristic (ROC) curve
}
\usage{
ROC(index, boolean, mask=NULL, nthres=NULL, thres=NULL, NAval=0, P=NA, Q=NA, 
progress=FALSE)
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
an optional integer indicating the number of equal-interval thresholds to be evaluated for the ROC curve. See Details below
}
  \item{thres}{
an optional numeric vector of thresholds to be evaluated for the ROC curve. See Details below
}
  \item{NAval}{
value for nodata (NA values) in the mask map
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
}
\details{
thresholds are calculated as the unique values of the index map after masking out NA values (default option), if neither nthres nor thres is provided. The default option can be time-consuming if the amount of unique values in the index map (after masking out NA values) is large (e.g., > 1000). In the latter case, the user may prefer to enter specified thresholds (with the \code{thres} argument), or to indicate the number of equal-interval thresholds to be evaluated for the ROC curve (with the \code{nthres} argument)
}
\value{
an object of class \code{Roc} containing the ROC table, the area under the curve (AUC), maximum AUC and minimum AUC
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
rocd <- ROC(index, boolean, mask)
rocd


# thresholds can also be defined by indicating the number of equal-interval thresholds 
rocd <- ROC(index, boolean, mask, nthres=10)
rocd

# A vector of thresholds can be used to define the thresholds
rocd <- ROC(index, boolean, mask, thres=seq(0, 100, by=10))
rocd
}
\dontrun{
index <- raster("Prob_Map2.rst")
boolean <- raster("Change_Map2b.rst")
mask <- raster("MASK3.rst")

# thresholds can also be defined by indicating the number of equal-interval thresholds 
rocd <- ROC(index, boolean, mask, nthres=100)
rocd

# A vector of thresholds can be used to define the thresholds
rocd <- ROC(index, boolean, mask, thres=seq(0, 100, by=10))
rocd
}

}

\keyword{ spatial }
