\name{quantile.dtn}
\alias{dquantile}
\alias{pquantile}
\alias{qquantile}
\alias{rquantile}

\title{
Interpolated quantile distribution with exponential tails
}
\description{
\code{dquantile} gives a probability density function (pdf) by combining
step-interpolation of probability densities for specified
\code{tau}-quantiles (\code{quant}) with exponential lower/upper tails
(Quinonero-Candela, 2006). Point mass (e.g., as might occur when using
censored QRNN models) can be defined by setting \code{lower} to the left
censoring point. \code{pquantile} gives the cumulative distribution
function (cdf); the \code{\link{integrate}} function is used for values
outside the range of \code{quant}. The inverse cdf is given by
\code{qquantile}; the \code{\link{uniroot}} function is used for values
outside the range of \code{tau}. \code{rquantile} is used for
generating random variates.
}
\usage{
dquantile(x, tau, quant, lower = -Inf)
pquantile(q, tau, quant, lower = -Inf, ...)
qquantile(p, tau, quant, lower = -Inf,
          tol = .Machine$double.eps^0.25, maxiter = 1000,
          range.mult = 10, max.error = 10, ...)
rquantile(n, tau, quant, lower = -Inf,
          tol = .Machine$double.eps^0.25, maxiter = 1000,
          range.mult = 10, max.error = 10, ...)
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
      additional arguments passed to \code{\link{integrate}}.
  }
}
\value{
\code{dquantile} gives the pdf, \code{pquantile} gives the cdf,
\code{qquantile} gives the inverse cdf (or quantile function), and
\code{rquantile} generates random deviates.
}
\references{
Quinonero-Candela, J., C. Rasmussen, F. Sinz, O. Bousquet,
B. Scholkopf, 2006. Evaluating Predictive Uncertainty Challenge.
Lecture Notes in Artificial Intelligence, 3944: 1-27.
}
\seealso{
\code{\link{integrate}}, \code{\link{uniroot}}, \code{\link{qrnn.predict}}
}
\examples{
## Normal distribution

tau <- seq(0.05, 0.95, by = 0.05)
quant <- qnorm(tau)

x <- seq(-3, 3, length = 500)
plot(x, dnorm(x), type = "l", col = "red", lty = 2, lwd = 2,
     main = "pdf")
lines(x, dquantile(x, tau, quant), col = "blue")

q <- seq(-3, 3, length = 20)
plot(q, pnorm(q), type = "b", col = "red", lty = 2, lwd = 2,
     main = "cdf")
lines(q, pquantile(q, tau, quant),
      col = "blue")

p <- c(0.001, 0.01, 0.025, seq(0.05, 0.95, by=0.05),
       0.975, 0.99, 0.999)
plot(p, qnorm(p), type = "b", col = "red", lty = 2, lwd = 2,
     main = "inverse cdf")
lines(p, qquantile(p, tau, quant), col = "blue")

## Distribution with point mass at zero

tau.0 <- c(0.3, 0.5, 0.7, 0.8, 0.9)
quant.0 <- c(0, 5, 7, 15, 20)

r.0 <- rquantile(500, tau = tau.0, quant = quant.0, lower = 0)
x.0 <- seq(0, 40, by = 0.5)
d.0 <- dquantile(x.0, tau = tau.0, quant = quant.0, lower = 0)
p.0 <- pquantile(x.0, tau = tau.0, quant = quant.0, lower = 0)
q.0 <- qquantile(p.0, tau = tau.0, quant = quant.0, lower = 0)

par(mfrow = c(2, 2))
plot(r.0, pch = 20, main = "random")
plot(x.0, d.0, type = "b", col = "red", main = "pdf")
plot(x.0, p.0, type = "b", col = "blue", ylim = c(0, 1),
     main = "cdf")
plot(p.0, q.0, type = "b", col = "green", xlim = c(0, 1),
     main = "inverse cdf")
}
