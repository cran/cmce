\name{gp.realize}
\alias{gp.realize}
\title{Draw an unconditional GP realization.}
\usage{
gp.realize(l.v,mu,lambda,rhos,lambdaf=Inf,eps=1e-10,from="")
}
\arguments{
\item{l.v}{A list of difference matrices for the design as calculated by \code{makedistlist()}}
\item{mu}{An \code{n x 1} mean vector for the realization}
\item{lambda}{A scalar quantity for the marginal precision of the drawn realization}
\item{rhos}{A \code{p x 1} vector of correlation parameters}
\item{lambdaf}{A scalar quantity denoting the precision of the error (i.i.d.) component of the realization}
\item{eps}{A fudge factor to help with inverting large correlation matrices}
\item{from}{Internal use only}
}
\value{
An \code{n x 1} vector of the GP realization calculated over the finite locations of the original design matrix.
}
\description{
The list \code{l.v} is created by \code{makedistlist()} from a {n x p} design matrix on which to generate the unconditional GP realization.  The correlation parameters take on values in \code{(0,1)} and use the Gaussian correlation model, calculated as \code{rho^d(xi,xj)^2}.  This helper function is used to create the demonstration data used in the example for \code{calibrate()}.
}
\examples{
library(cmce)

design=matrix(runif(10,1,5),ncol=2,nrow=5)
r=getranges(design)
design=scaledesign(design,r)
l.v=makedistlist(design)
rho=c(0.2,0.01)
muv=rep(0,nrow(design))
lambdav=1
surface=gp.realize(l.v,muv,lambdav,rho)

}
\seealso{
  \code{\link{getranges}}
  \code{\link{scaledesign}}
  \code{\link{makedistlist}}
}
