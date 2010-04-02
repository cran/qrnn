\name{hramp.prime}
\alias{hramp.prime}
\title{
Derivative of the Huber norm approximation to the ramp function
}
\description{
Derivative of the Huber norm approximation to the ramp function.
}
\usage{
hramp.prime(x, lower, eps)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{lower}{
      left censoring point.
  }
  \item{eps}{
      epsilon value used in \code{\link{huber}} and related functions.
  }  
}
\seealso{
\code{\link{hramp}}, \code{\link{huber}}
}
\examples{
x <- seq(-10, 10, length = 100)
plot(x, hramp.prime(x, lower = 0, eps = 1))
}
