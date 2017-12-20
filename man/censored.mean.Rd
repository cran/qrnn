\name{censored.mean}
\alias{censored.mean}
\title{
A hybrid mean/median function for left censored variables
}
\description{
Returns the median if the majority of values are censored and the mean otherwise.
}
\usage{
censored.mean(x, lower, trim=0)
}
\arguments{
  \item{x}{
      numeric vector.
  }
  \item{lower}{
      left censoring point.
  }
  \item{trim}{
       fraction of observations to be trimmed from each end of \code{x} before the mean is computed.
  }
}
\seealso{
\code{\link{qrnn.fit}}, \code{\link{qrnn.predict}}}
\examples{
x <- c(0, 0, 1, 2, 3)
print(censored.mean(x, lower=0))
x.cens <- c(0, 0, 0, 1, 2)
print(censored.mean(x.cens, lower=0))
}
