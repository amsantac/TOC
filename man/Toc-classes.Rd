\name{Toc-class}

\docType{class}

\alias{Toc-class}
\alias{Roc-class}

\title{Toc and Roc classes}

\description{
Toc and Roc classes
}


\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Toc", ...)}, or with the helper functions such as \code{Toc}.
}

\section{Slots}{
Slots for Roc and Toc objects
	\describe{
    \item{\code{table}:}{data.frame} 
    \item{\code{AUC}:}{numeric; Area Under the Curve}
    \item{\code{maxAUC}:}{numeric; maximum AUC}
    \item{\code{minAUC}:}{numeric; minimum AUC}
    \item{\code{prevalence}:}{numeric; prevalence}
    \item{\code{population}:}{numeric; population}
    \item{\code{units}:}{character; units for data in the TOC table, prevalence and population} 
  }
}

\examples{
showClass("Toc")
}

\keyword{classes}
\keyword{spatial}
