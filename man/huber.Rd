\name{huber}
\alias{huber}
\title{
Huber norm function
}
\description{
Huber norm function providing a hybrid L1/L2 norm
}
\usage{
huber(x, eps)
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
\code{\link{huber.prime}}
}
\examples{
x <- seq(-10, 10, length = 100)
plot(x, huber(x, eps = 1))
}
