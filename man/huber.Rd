\name{huber}
\alias{huber}
\alias{huber.prime}
\alias{hramp}
\alias{hramp.prime}
\alias{tilted.approx}
\alias{tilted.approx.prime}
\title{
Huber norm and Huber approximations to the ramp and tilted absolute value functions
}
\description{
Huber norm function providing a hybrid L1/L2 norm. Huber approximations to the ramp \code{hramp} and tilted absolute value \code{tilted.approx} functions. \code{huber.prime}, \code{hramp.prime}, and \code{tilted.approx.prime} provide the corresponding derivatives.
}
\usage{
huber(x, eps)
huber.prime(x, eps)
hramp(x, lower, eps)
hramp.prime(x, lower, eps)
tilted.approx(x, tau, eps)
tilted.approx.prime(x, tau, eps)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{eps}{
      epsilon value used in \code{\link{huber}} and related functions.
  }
  \item{tau}{
      desired tau-quantile.
  }
  \item{lower}{
      left censoring point.
  }
}
\seealso{
\code{\link{tilted.abs}}, \code{\link{qrnn.cost}}
}
\examples{
x <- seq(-10, 10, length=100)
plot(x, huber(x, eps=1), type="l", col="black", ylim=c(-2, 10), ylab="")
lines(x, hramp(x, lower=0, eps=1), col="red")
lines(x, tilted.approx(x, tau=0.1, eps=1), col="blue")
lines(x, huber.prime(x, eps=1), col="black", lty=2)
lines(x, hramp.prime(x, lower=0, eps=1), lty=2, col="red")
lines(x, tilted.approx.prime(x, tau=0.1, eps=1), lty=2, col="blue")
}
