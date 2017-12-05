\name{softplus.prime}
\alias{softplus.prime}
\title{
Derivative of the approximation to the rectified linear function
}
\description{
Derivative of the approximation to the rectified linear function.
}
\usage{
softplus.prime(x, alpha=2)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{alpha}{
      rounded transition parameter.
  }
}
\seealso{
\code{\link{softplus}}
}
\examples{
x <- seq(-2, 2, length=200)
plot(x, softplus.prime(x, alpha=2))
}
