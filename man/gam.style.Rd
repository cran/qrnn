\name{gam.style}
\alias{gam.style}
\title{
GAM-style effects plots for interpreting qrnn models
}
\description{
GAM-style effects plots provide a graphical means of interpreting
relationships between predictors and conditional quantiles
predicted by a qrnn. From Plate et al. (2000): The effect of the
\code{i}th input variable at a particular input point \code{Delta.i.x}
is the change in \code{f} resulting from changing \code{X1} to \code{x1}
from \code{b1} (the baseline value [...]) while keeping the other
inputs constant. The effects are plotted as short line segments, centered
at (\code{x.i}, \code{Delta.i.x}), where the slope of the segment
is given by the partial derivative. Variables that strongly influence
the function value have a large total vertical range of effects.
Functions without interactions appear as possibly broken straight lines
(linear functions) or curves (nonlinear functions). Interactions show up as
vertical spread at a particular horizontal location, that is, a vertical
scattering of segments. Interactions are present when the effect of
a variable depends on the values of other variables.
}

\usage{
gam.style(x, parms, column, baseline = mean(x[,column]),
         epsilon = 1e-5, seg.len = 0.02, seg.cols = "black",
         plot = TRUE, return.results = FALSE, trim = 0,
         ...)
}
\arguments{
  \item{x}{
      matrix with number of rows equal to the number of samples and number of columns equal to the number of predictor variables.
}
  \item{parms}{
      list returned by \code{\link{qrnn.fit}}.
}
  \item{column}{
      column of \code{x} for which effects plots should be returned.
}
  \item{baseline}{
      value of \code{x[,column]} to be used as the baseline for calculation of predictor effects; defaults to \code{mean(x[,column])}.
}
  \item{epsilon}{
      step-size used in the finite difference calculation of the partial derivatives.
}
  \item{seg.len}{
      length of effects line segments expressed as a fraction of the range of \code{x[,column]}.
}
  \item{seg.cols}{
      colors of effects line segments.
}
  \item{plot}{
      if \code{TRUE} (the default) then an effects plots for the given model is produced.
}
  \item{return.results}{
      if \code{TRUE} then values of effects and partial derivatives are returned.
}
  \item{trim}{
     if \code{plot=TRUE} and \code{parms} is for a model with \code{n.ensemble > 1}, value of \code{trim} passed to \code{\link{censored.mean}}.
}
  \item{\dots}{
      further arguments to be passed to \code{plot}.
}
}
\seealso{
    \code{\link{qrnn.fit}}, \code{\link{qrnn.predict}}
}
\value{
A list with elements:
  \item{effects}{a matrix of predictor effects.}
  \item{partials}{a matrix of predictor partial derivatives.}
}
\examples{
## YVR precipitation data with seasonal cycle and NCEP/NCAR Reanalysis
## covariates

data(YVRprecip)

y <- YVRprecip$precip
x <- cbind(sin(2*pi*seq_along(y)/365.25),
           cos(2*pi*seq_along(y)/365.25),
           YVRprecip$ncep)

## Fit QRNN and quantile regression models for the 95th percentile

train <- c(TRUE, rep(FALSE, 9))
w.qrnn <- qrnn.fit(x = x[train,], y = y[train,,drop=FALSE],
                   n.hidden = 2, tau = 0.95, iter.max = 1000,
                   n.trials = 1, lower = 0)
w.qreg <- qrnn.fit(x = x[train,], y = y[train,,drop=FALSE],
                   tau = 0.95, iter.max = 1000, n.trials = 1,
                   lower = 0, Th = linear, Th.prime = linear.prime)

## GAM-style plots for slp, sh700, and z500

for (column in 3:5) {
    gam.style(x[train,], parms = w.qrnn, column = column,
              main = "QRNN")
    gam.style(x[train,], parms = w.qreg, column = column,
              main = "QREG")
}
}
\references{
    Cannon, A.J. and I.G. McKendry, 2002. A graphical sensitivity analysis
    for interpreting statistical climate models: Application to Indian
    monsoon rainfall prediction by artificial neural networks and
    multiple linear regression models. International Journal of
    Climatology, 22:1687-1708.
    
    Plate, T., J. Bert, J. Grace, and P. Band, 2000. Visualizing the function
    computed by a feedforward neural network. Neural Computation,
    12(6): 1337-1354.
}
