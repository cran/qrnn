\name{hramp.prime}
\alias{hramp.prime}
\title{
Derivative of the approximation to the ramp function
}
\description{
Derivative of the approximation to the ramp function.
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
      epsilon value used in the approximation.
  }  
}
\seealso{
\code{\link{hramp}}, \code{\link{huber}}
}
\examples{
x <- seq(-2, 2, length=200)
plot(x, hramp.prime(x, lower=0, eps=1))
}
