\name{tilted.huber}
\alias{tilted.huber}
\title{
Tilted Huber norm function
}
\description{
Huber norm approximation to the tilted absolute value function.
}
\usage{
tilted.huber(x, tau, eps)
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
\code{\link{huber}}, \code{\link{tilted.abs}},
\code{\link{tilted.huber.prime}}
}
\examples{
x <- seq(-10, 10, length = 100)
plot(x, tilted.huber(x, tau = 0.75, eps = 1))
}
