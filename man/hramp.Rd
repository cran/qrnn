\name{hramp}
\alias{hramp}
\title{
Huber norm approximation to the ramp function
}
\description{
Huber norm approximation to the ramp function.
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
      epsilon value used in \code{\link{huber}} and related functions.
  }  
}
\seealso{
\code{\link{hramp.prime}}, \code{\link{huber}}
}
\examples{
x <- seq(-10, 10, length = 100)
plot(x, hramp(x, lower = 0, eps = 1))
}
