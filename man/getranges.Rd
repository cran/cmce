\name{getranges}
\alias{getranges}
\title{Get variable ranges from a design matrix.}
\usage{
getranges(design)
}
\arguments{
\item{design}{An \code{n x p} matrix of input settings}

}
\value{
A \code{p x 2} matrix with the lower and upper bounds (rounded to nearest integer value) of all \code{p} variables in the design matrix.
}
\description{
\code{getranges()} is a helper function to get the lower/upper bounds of variables in a design matrix, used for rescaling the inputs to the \code{[0,1])} hypercube.
}
\examples{
library(cmce)

design=matrix(runif(10,1,5),ncol=2,nrow=5)
getranges(design)
}
