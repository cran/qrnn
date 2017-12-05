\name{softplus}
\alias{softplus}
\title{
Approximation to the rectified linear function
}
\description{
Approximation to the rectified linear function. Used as the hidden layer transfer function for nonlinear QRNN models.
}
\usage{
softplus(x, alpha=2)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{alpha}{
      rounded transition parameter.
  }
}
\seealso{
\code{\link{softplus.prime}}
}
\examples{
x <- seq(-2, 2, length=200)
plot(x, softplus(x, alpha=2))
}
