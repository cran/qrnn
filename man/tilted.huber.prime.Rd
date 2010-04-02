\name{tilted.huber.prime}
\alias{tilted.huber.prime}
\title{
Derivative of the tilted Huber norm function
}
\description{
Derivative of the Huber norm approximation to the tilted absolute value function.
}
\usage{
tilted.huber.prime(x, tau, eps)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{tau}{
      desired tau-quantile.
  }
  \item{eps}{
      epsilon value used in \code{\link{huber}} and related functions.
  }  
}
\seealso{
\code{\link{tilted.huber}}, \code{\link{huber}}, \code{\link{tilted.abs}}
}
\examples{
x <- seq(-10, 10, length = 100)
plot(x, tilted.huber.prime(x, tau = 0.75, eps = 1))
}
