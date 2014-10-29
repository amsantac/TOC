\name{TOC}
\alias{TOC}

\title{
calculates TOC table
}
\description{
calculates TOC
}
\usage{
TOC(input, ref, mask=NULL, nthres, NAvalue=0)
}

\arguments{
  \item{input}{
input Raster map
}
  \item{ref}{
reference Raster map
}
  \item{mask}{
mask Raster map
}
  \item{nthres}{
integer. Number of thresholds
}
  \item{NAvalue}{
value for nodata
}
}
\value{
a data.frame with TOC values
}
\seealso{
\code{\link{TOCplot}}
}
\examples{
input <- raster(system.file("external/Prob_Map2.rst", package="composite"))
ref <- raster(system.file("external/Change_Map2b.rst", package="composite"))
mask <- raster(system.file("external/MASK3.rst", package="composite"))
tocd <- TOC(input, ref, mask, nthres=100, NAvalue=0)
}

\keyword{ spatial }
