\name{linear}
\alias{linear}
\title{
Identity function
}
\description{
Computes a trivial identity function. Used as the hidden layer transfer function for linear QRNN models.
}
\usage{
linear(x)
}
\arguments{
  \item{x}{
      numeric vector.
  }
}
\seealso{
\code{\link{linear.prime}}, \code{\link{sigmoid}}
}
\examples{
x <- seq(-10, 10, length=100)
plot(x, linear(x))
}
