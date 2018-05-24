\name{cmce-package}
\alias{cmce-package}
\docType{package}
\title{
  \code{cmce} is an open-source R package implementing Bayesian calibration for deterministic and stochastic simulation models.
}
\description{
  \code{cmce} implements a Bayesian calibration model based on a dimension-reduction approach via empirical orthogonal functions.  The model generalizes the popular SVD-based approach (e.g. Higdon et al. (2008)) to a tensor-variate decomposition, enabling the calibration of both deterministic simulators and stochastic simulators.  Currently, \code{cmce} is a proof of concept implementation based entirely on \code{R} code, leaving open the possibility of future performance improvements if the codebase is moved to \code{C} at a later point.  It supports both the popular additive discrepancy as well as multiplicative discrepancy.  The model is fit using a Metropolis-within-Gibbs MCMC algorithm.
}
\details{
  The main model fitting function in the package is \code{calibrate()}.  This function will return posterior samples of the calibration parameters as well as the predicted field process, predicted discrepancies and predicted calibrated simulator.  The model currently does not support predicting off the input grid and assumes that field data and simulator outputs are observed on the same input grid.  
}
\author{
Matthew T. Pratola <mpratola@stat.osu.edu> [aut, cre, cph]
}
\references{
  Pratola, Matthew T. and Chrebtii, Oksana. (2018)
  Bayesian Calibration of Multistate Stochastic Simulators. 
  \emph{Statistica Sinica}, \bold{28}, 693--720.
  \doi{10.5705/ss.202016.0403}.

  Higdon, Dave, Gattiker, James, Williams, Brian and Rightley, Maria. (2008)
  Computer model calibration using high-dimensional output. 
  \emph{Journal of the American Statistical Association}, \bold{103}, 570--583.
  \doi{10.1198/016214507000000888}.
}
\keyword{ package }
\seealso{
  \code{\link{calibrate}}
}
