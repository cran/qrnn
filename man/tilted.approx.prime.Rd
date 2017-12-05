\name{tilted.approx.prime}
\alias{tilted.approx.prime}
\title{
Derivative of the approximate tilted absolute value function
}
\description{
Derivative of the approximate to tilted absolute value function.
}
\usage{
tilted.approx.prime(x, tau, eps)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{tau}{
      desired tau-quantile.
  }
  \item{eps}{
      epsilon value used in the approximation.
  }  
}
\seealso{
\code{\link{tilted.approx}}, \code{\link{tilted.abs}},
\code{\link{huber}}
}
\examples{
x <- seq(-2, 2, length=200)
plot(x, tilted.approx.prime(x, tau=0.75, eps=1))
}
