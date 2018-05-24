\name{unscalemat}
\alias{unscalemat}
\title{Unscale a matrix back to its original ranges.}
\usage{
unscalemat(mat,r)
}
\arguments{
\item{mat}{An \code{n x p} matrix of numbers scaled to the \code{[0,1]} hypercube}
\item{r}{An \code{p x 2} matrix of the original ranges of the variables}

}
\value{
A \code{n x p} matrix with variables rescaled back to their original ranges, as specified by \code{ranges}.
}
\description{
\code{unscalemat()} is a helper function to rescale a matrix back to its original ranges.  Typically this is used to rescale the posterior samples of the parameters back to their original scale.
}
\examples{
library(cmce)

design=matrix(runif(10,1,5),ncol=2,nrow=5)
r=getranges(design)
design=scaledesign(design,r)
unscalemat(design,r)

}
\seealso{
  \code{\link{getranges}}
  \code{\link{scaledesign}}
}
