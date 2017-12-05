\name{tilted.approx}
\alias{tilted.approx}
\title{
Approximation to the tilted absolute value function
}
\description{
Approximation to the tilted absolute value function.
}
\usage{
tilted.approx(x, tau, eps)
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
\code{\link{tilted.abs}}, \code{\link{tilted.approx.prime}},
\code{\link{huber}}
}
\examples{
x <- seq(-2, 2, length=200)
plot(x, tilted.approx(x, tau=0.75, eps=1))
}
