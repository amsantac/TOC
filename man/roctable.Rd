\name{roctable}
\alias{roctable}

\title{
Construct a basic ROC table
}
\description{
TOC internal function. Construct a basic ROC table
}
\usage{
roctable(indval, boolval, maskval = NULL, nthres = NULL, 
         thres = NULL, NAval = 0, progress = FALSE, 
         ones.bool = NULL, zeros.bool = NULL)
}

\arguments{
  \item{indval}{
numeric index vector
}
  \item{boolval}{
numeric boolean vector
}
  \item{maskval}{
numeric mask vector
}
  \item{nthres}{
an optional integer indicating the number of equal-interval thresholds to be evaluated for the TOC curve. See Details below
}
  \item{thres}{
an optional numeric vector of thresholds to be evaluated for the TOC curve. See Details below
}
  \item{NAval}{
value for nodata (NA values) in the mask map
}
  \item{progress}{
logical; if TRUE, a progress bar is shown
}
  \item{ones.bool}{
numeric value indicating total number of 1's in the boolean vector
}
  \item{zeros.bool}{
numeric value indicating total number of 0's in the boolean vector
}
}
\note{
This function is not meant to be called by users directly
}
\value{
a data.frame with a basic ROC table and a numeric value for minimum value in the index vector
}

\keyword{ spatial }
