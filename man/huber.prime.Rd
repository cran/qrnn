\name{huber.prime}
\alias{huber.prime}
\title{
Derivative of the Huber norm function
}
\description{
Derivative of the Huber norm function.
}
\usage{
huber.prime(x, eps)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{eps}{
      epsilon value used in \code{\link{huber}} and related functions.
  }  
}
\seealso{
\code{\link{huber}}
}
\examples{
x <- seq(-10, 10, length = 100)
plot(x, huber.prime(x, eps = 1))
}
