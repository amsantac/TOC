\name{TOC-package}
\alias{TOC-package}
\docType{package}
\title{
Total Operating Characteristic (TOC) Curve and ROC Curve
}
\description{
Construction of the Total Operating Characteristic (TOC) Curve and the Receiver (aka Relative) Operating Characteristic (ROC) Curve for spatial and non-spatial data. The TOC method is a modification of the ROC method which measures the ability of an index variable to diagnose either presence or absence of a characteristic. The diagnosis depends on whether the value of an index variable is above a threshold. Each threshold generates a two-by-two contingency table, which contains four entries: hits (H), misses (M), false alarms (FA), and correct rejections (CR). While ROC shows for each threshold only two ratios, H/(H + M) and FA/(FA + CR), TOC reveals the size of every entry in the contingency table for each threshold (Pontius Jr., R.G., Si, K. 2014. The total operating characteristic to measure diagnostic ability for multiple thresholds. Int. J. Geogr. Inf. Sci. 28 (3), 570-583 <doi:10.1080/13658816.2013.862623>).}
\details{
\tabular{ll}{
Package: \tab TOC\cr
Type: \tab Package\cr
Version: \tab 0.0-5\cr
Date: \tab 2020-05-04\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
}
\author{
Robert G. Pontius, Ali Santacruz, Amin Tayyebi, Benoit Parmentier, Kangping Si

Maintainer: Ali Santacruz
}
\references{
Pontius Jr., R.G., Kangpin, Si. 2014. \emph{The total operating characteristic to measure diagnostic ability for multiple thresholds}. International Journal of Geographical Information Science 28 (3): 570-583. <doi:10.1080/13658816.2013.862623>

Pontius, G., Parmentier, B. 2014. \emph{Recommendations for using the Relative Operating Characteristic (ROC)}. Landscape Ecology 29 (3): 367-382. <doi:10.1007/s10980-013-9984-8>
}
\seealso{
\code{\link{TOC}}, \code{\link{plot}} 
}
\keyword{ package }
\keyword{ spatial }