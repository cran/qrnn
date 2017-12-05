\name{elu.prime}
\alias{elu.prime}
\title{
Derivative of the exponential linear function
}
\description{
Derivative of the exponential linear function.
}
\usage{
elu.prime(x, alpha=1)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{alpha}{
      transition parameter.
  }
}
\seealso{
\code{\link{elu}}
}
\examples{
x <- seq(-2, 2, length=200)
plot(x, elu.prime(x, alpha=1))
}
