\name{elu}
\alias{elu}
\title{
Exponential linear function
}
\description{
Exponential linear function. Used as the hidden layer transfer function for nonlinear QRNN models.
}
\usage{
elu(x, alpha=1)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{alpha}{
      transition parameter.
  }
}
\seealso{
\code{\link{elu.prime}}
}
\examples{
x <- seq(-2, 2, length=200)
plot(x, elu(x, alpha=1))
}
