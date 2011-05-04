\name{qrnn-package}
\alias{qrnn-package}
\alias{qrnn}
\docType{package}
\title{Quantile regression neural network}
\description{
This package implements the quantile regression neural network (QRNN)
(Taylor, 2000), which is the artificial neural network  analog of
linear quantile regression. The QRNN formulation follows from previous
work on the estimation of censored regression quantiles, thus allowing
predictions for mixed discrete-continuous variables like precipitation
(Friederichs and Hense, 2007). A differentiable
approximation to the quantile regression cost function is adopted so
that a simplified form of the finite smoothing algorithm (Chen, 2007) can be
used to estimate model parameters.

QRNN models are fit using the \code{\link{qrnn.fit}} function. Depending
on the choice of hidden layer transfer function, one can also fit
quantile (ridge) regression models. A kernel quantile ridge regression
model can be specified with the aid of the \code{\link{qrnn.rbf}}
function.

Predictions from a fitted model are made using the
\code{\link{qrnn.predict}} function. If models for multiple quantiles
have been fitted, \code{\link{dquantile}} and its companion functions
are available to approximate a probability density function and related
distribution functions.

With the exception of \code{\link{gam.style}}, which can be used
to investigate fitted predictor/predictand relationships, most other functions
are used internally and should not need to be called directly by the user.
}
\details{
\tabular{ll}{
Package: \tab qrnn\cr
Type: \tab Package\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Alex J. Cannon <acannon@eos.ubc.ca>
}
\references{

Cannon, A.J., 2010. Quantile regression neural networks: implementation
in R and application to precipitation downscaling. Computers & Geosciences,
doi:10.1016/j.cageo.2010.07.005

Chen, C., 2007. A finite smoothing algorithm for quantile regression.
Journal of Computational and Graphical Statistics, 16: 136-164.

Friederichs, P. and A. Hense, 2007. Statistical downscaling of extreme
precipitation events using censored quantile regression. Monthly Weather
Review, 135: 2365-2378. 

Quinonero-Candela, J., C. Rasmussen, F. Sinz, O. Bousquet,
B. Scholkopf, 2006. Evaluating Predictive Uncertainty Challenge.
Lecture Notes in Artificial Intelligence, 3944: 1-27.

Taylor, J.W., 2000.  A quantile regression neural network approach to
estimating the conditional density of multiperiod returns. Journal of
Forecasting, 19(4): 299-311.
}
