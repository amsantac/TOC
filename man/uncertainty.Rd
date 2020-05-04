\name{uncertainty}
\alias{uncertainty}

\title{
Uncertainty in AUC calculation
}
\description{
TOC internal function. It calculates uncertainty in AUC calculation
}
\usage{
uncertainty(index, tocd)
}

\arguments{
  \item{index}{
index vector
}
  \item{tocd}{
data.frame output from \code{roctable}
}
}
\note{
This function is not meant to be called by users directly
}
\value{
a numeric value representing uncertainty in AUC calculation 
}

\keyword{ spatial }
