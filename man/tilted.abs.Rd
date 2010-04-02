\name{tilted.abs}
\alias{tilted.abs}
\title{
Tilted absolute value function
}
\description{
Tilted absolute value function. Also known as the check function, hinge function, or the pinball loss function.
}
\usage{
tilted.abs(x, tau)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{tau}{
      desired tau-quantile.
  }
}
\seealso{
\code{\link{huber}}, \code{\link{tilted.huber}}}
\examples{
x <- seq(-10, 10, length = 100)
plot(x, tilted.abs(x, tau = 0.75))
}
