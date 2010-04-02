\name{sigmoid.prime}
\alias{sigmoid.prime}
\title{
Derivative of the sigmoid function
}
\description{
Derivative of the hyperbolic tangent sigmoid function.
}
\usage{
sigmoid.prime(x)
}
\arguments{
  \item{x}{
      numeric vector.
  }  
}
\seealso{
\code{\link{sigmoid}}, \code{\link{linear}}
}
\examples{
x <- seq(-10, 10, length = 100)
plot(x, sigmoid.prime(x))
}
