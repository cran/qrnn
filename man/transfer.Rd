\name{transfer}
\alias{sigmoid}
\alias{sigmoid.prime}
\alias{elu}
\alias{elu.prime}
\alias{softplus}
\alias{softplus.prime}
\alias{linear}
\alias{linear.prime}
\title{
Transfer functions and their derivatives
}
\description{
The \code{sigmoid}, exponential linear \code{elu}, and
\code{softplus} functions can be used as the hidden layer
transfer function for a nonlinear QRNN model. \code{sigmoid} is
used by default. The \code{linear} function is used as the
hidden layer transfer function for linear QRNN models.
\code{sigmoid.prime}, \code{elu.prime}, \code{softplus.prime},
and \code{linear.prime} provide the corresponding derivatives.
}
\usage{
sigmoid(x)
sigmoid.prime(x)
elu(x, alpha=1)
elu.prime(x, alpha=1)
softplus(x, alpha=2)
softplus.prime(x, alpha=2)
linear(x)
linear.prime(x)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{alpha}{
      transition parameter for \code{elu} and \code{softplus} functions.
  }
}
\examples{
x <- seq(-10, 10, length=100)
plot(x, sigmoid(x), type="l", col="black", ylab="")
lines(x, sigmoid.prime(x), lty=2, col="black")
lines(x, elu(x), col="red")
lines(x, elu.prime(x), lty=2, col="red")
lines(x, softplus(x), col="blue")
lines(x, softplus.prime(x), lty=2, col="blue")
lines(x, linear(x), col="green")
lines(x, linear.prime(x), lty=2, col="green")
}
