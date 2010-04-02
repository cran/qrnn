\name{linear.prime}
\alias{linear.prime}
\title{
Derivitative of the linear function
}
\description{
Derivative of the linear function.
}
\usage{
linear.prime(x)
}
\arguments{
  \item{x}{
      numeric vector.
  }
}
\seealso{
\code{\link{linear}}, \code{\link{sigmoid}}
}
\examples{
x <- seq(-10, 10, length = 100)
plot(x, linear.prime(x))
}
