\name{hramp}
\alias{hramp}
\title{
Approximation to the ramp function
}
\description{
Approximation to the ramp function.
}
\usage{
hramp(x, lower, eps)
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
\code{\link{hramp.prime}}, \code{\link{huber}}
}
\examples{
x <- seq(-2, 2, length=200)
plot(x, hramp(x, lower=0, eps=1))
}
