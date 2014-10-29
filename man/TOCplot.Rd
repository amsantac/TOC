\name{TOCplot}
\alias{TOCplot}

\title{
generates TOC plot
}
\description{
generates TOC plot
}
\usage{
TOCplot(truePositive, truePositRate, falsePositive, falsePositRate, population)
}

\arguments{
  \item{truePositive}{
true positives
}
  \item{falsePositive}{
false positives
}
  \item{truePositRate}{
true positive rates
}
  \item{falsePositRate}{
false positive rates
}
  \item{population}{
population
}
}
\value{
a TOC plot
}
\seealso{
\code{\link{TOC}}
}
\examples{
input <- raster(system.file("external/Prob_Map2.rst", package="composite"))
ref <- raster(system.file("external/Change_Map2b.rst", package="composite"))
mask <- raster(system.file("external/MASK3.rst", package="composite"))
tocd <- TOC(input, ref, mask, nthres=100, NAvalue=0)

ncells <- 1270089
cellSize <- 1000
population <- (ncells * cellSize^2)/1000^2
TOCplot(tocd$truePositive, tocd$truePosRate, tocd$falsePositive, tocd$falsePosRate, population)
}

\keyword{ spatial }
