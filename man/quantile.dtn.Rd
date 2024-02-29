\name{quantile.dtn}
\alias{dquantile}
\alias{pquantile}
\alias{qquantile}
\alias{rquantile}
\alias{pquantile.nw}
\alias{qquantile.nw}
\alias{rquantile.nw}
\encoding{UTF-8}
\title{
Interpolated quantile distribution with exponential tails and Nadaraya-Watson quantile distribution
}
\description{
\code{dquantile} gives a probability density function (pdf) by combining
step-interpolation of probability densities for specified
\code{tau}-quantiles (\code{quant}) with exponential lower/upper tails
(Quiñonero-Candela, 2006; Cannon, 2011). Point mass (e.g., as might occur
when using censored QRNN models) can be defined by setting \code{lower} to
the left censoring point. \code{pquantile} gives the cumulative distribution
function (cdf); the \code{\link{integrate}} function is used for values
outside the range of \code{quant}. The inverse cdf is given by
\code{qquantile}; the \code{\link{uniroot}} function is used for values
outside the range of \code{tau}. \code{rquantile} is used for
generating random variates. 

Alternative formulations (without left censoring) based on the
Nadaraya-Watson estimator \code{[p,q,r]quantile.nw} are also provided
(Passow and Donner, 2020).

Note: these functions have not been extensively tested or optimized and 
should be considered experimental.

}
\usage{
dquantile(x, tau, quant, lower=-Inf)
pquantile(q, tau, quant, lower=-Inf, ...)
pquantile.nw(q, tau, quant, h=0.001, ...)
qquantile(p, tau, quant, lower=-Inf,
          tol=.Machine$double.eps^0.25, maxiter=1000,
          range.mult=1.1, max.error=100, ...)
qquantile.nw(p, tau, quant, h=0.001)
rquantile(n, tau, quant, lower=-Inf,
          tol=.Machine$double.eps^0.25, maxiter=1000,
          range.mult=1.1, max.error=100, ...)
rquantile.nw(n, tau, quant, h=0.001)
}
\arguments{
  \item{x, q}{
      vector of quantiles.
  }
  \item{p}{
      vector of cumulative probabilities.
  }
  \item{n}{
      number of random samples.
  }
  \item{tau}{
      ordered vector of cumulative probabilities associated with \code{quant} argument.
  }
  \item{quant}{
      ordered vector of quantiles associated with \code{tau} argument.
  }
  \item{lower}{
      left censoring point.
  }
  \item{tol}{
      tolerance passed to \code{\link{uniroot}}.
  }
  \item{h}{
      bandwidth for Nadaraya-Watson kernel.
  }
  
  \item{maxiter}{
      maximum number of iterations passed to \code{\link{uniroot}}.
  }
  \item{range.mult}{
      values of \code{lower} and \code{upper} in \code{\link{uniroot}} are initialized to \cr \code{quant[1]-range.mult*diff(range(quant))} and \cr \code{quant[length(quant)]+range.mult*diff(range(quant))} respectively; \code{range.mult} is squared, \code{lower} and \code{upper} are recalculated, and \code{\link{uniroot}} is rerun if the current values lead to an exception.
  }
  \item{max.error}{
      maximum number of \code{uniroot} errors allowed before termination.
  }
  \item{...}{
      additional arguments passed to \code{\link{integrate}} or \code{\link{uniroot}}.
  }
}
\value{
\code{dquantile} gives the pdf, \code{pquantile} gives the cdf,
\code{qquantile} gives the inverse cdf (or quantile function), and
\code{rquantile} generates random deviates.
}
\references{
Cannon, A.J., 2011. Quantile regression neural networks: implementation
in R and application to precipitation downscaling. Computers & Geosciences,
37: 1277-1284. doi:10.1016/j.cageo.2010.07.005

Passow, C., R.V. Donner, 2020. Regression-based distribution mapping for
bias correction of climate model outputs using linear quantile regression.
Stochastic Environmental Research and Risk Assessment, 34:87-102.
doi:10.1007/s00477-019-01750-7

Quiñonero-Candela, J., C. Rasmussen, F. Sinz, O. Bousquet,
B. Scholkopf, 2006. Evaluating Predictive Uncertainty Challenge.
Lecture Notes in Artificial Intelligence, 3944: 1-27.
}
\seealso{
\code{\link{integrate}}, \code{\link{uniroot}}, \code{\link{qrnn.predict}}
}
\examples{
## Normal distribution
tau <- c(0.01, seq(0.05, 0.95, by=0.05), 0.99)
quant <- qnorm(tau)

x <- seq(-3, 3, length=500)
plot(x, dnorm(x), type="l", col="red", lty=2, lwd=2,
     main="pdf")
lines(x, dquantile(x, tau, quant), col="blue")

q <- seq(-3, 3, length=20)
plot(q, pnorm(q), type="b", col="red", lty=2, lwd=2,
     main="cdf")
lines(q, pquantile(q, tau, quant),
      col="blue")

abline(v=1.96, lty=2)
abline(h=pnorm(1.96), lty=2)
abline(h=pquantile(1.96, tau, quant), lty=3)
abline(h=pquantile.nw(1.96, tau, quant, h=0.01), lty=3)

p <- c(0.001, 0.01, 0.025, seq(0.05, 0.95, by=0.05),
       0.975, 0.99, 0.999)
plot(p, qnorm(p), type="b", col="red", lty=2, lwd=2,
     main="inverse cdf")
lines(p, qquantile(p, tau, quant), col="blue")

## Distribution with point mass at zero
tau.0 <- c(0.3, 0.5, 0.7, 0.8, 0.9)
quant.0 <- c(0, 5, 7, 15, 20)

r.0 <- rquantile(500, tau=tau.0, quant=quant.0, lower=0)
x.0 <- seq(0, 40, by=0.5)
d.0 <- dquantile(x.0, tau=tau.0, quant=quant.0, lower=0)
p.0 <- pquantile(x.0, tau=tau.0, quant=quant.0, lower=0)
q.0 <- qquantile(p.0, tau=tau.0, quant=quant.0, lower=0)

par(mfrow=c(2, 2))
plot(r.0, pch=20, main="random")
plot(x.0, d.0, type="b", col="red", main="pdf")
plot(x.0, p.0, type="b", col="blue", ylim=c(0, 1),
     main="cdf")
plot(p.0, q.0, type="b", col="green", xlim=c(0, 1),
     main="inverse cdf")
}
