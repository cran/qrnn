\name{sigmoid}
\alias{sigmoid}
\title{
Sigmoid function
}
\description{
Computes the hyperbolic tangent sigmoid function. Used as the hidden layer transfer function for nonlinear QRNN models.
}
\usage{
sigmoid(x)
}
\arguments{
  \item{x}{
      numeric vector.
  }
}
\seealso{
\code{\link{sigmoid.prime}}, \code{\link{linear}}
}
\examples{
x <- seq(-10, 10, length = 100)
plot(x, sigmoid(x))
}
