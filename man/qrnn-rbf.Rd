\name{qrnn.rbf}
\alias{qrnn.rbf}
\title{
Radial basis function kernel
}
\description{
Evaluate a kernel matrix based on the radial basis function kernel. Can
be used in conjunction with \code{\link{qrnn.fit}} with \code{Th} set to
\code{\link{linear}} and \code{penalty} set to a nonzero value for
kernel quantile ridge regression.
}
\usage{
qrnn.rbf(x, x.basis, sigma)
}
\arguments{
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
  }
  \item{x.basis}{
      covariate matrix with number of rows equal to the number of basis functions and number of columns equal to the number of variables.
  }
  \item{sigma}{
      kernel width
  }
}
\value{
    kernel matrix with number of rows equal to the number of samples and number of columns equal to the number of basis functions.
}
\seealso{
    \code{\link{qrnn.fit}}
}
\examples{
data(sinc)
x <- sinc$x
y <- sinc$y

set.seed(1)
kern <- qrnn.rbf(x, x.basis = x, sigma = 1)

parms <- qrnn.fit(x = kern, y = y, tau = 0.5, penalty = 0.1,
                  Th = linear, Th.prime = linear.prime,
                  iter.max = 500, n.trials = 1)
p <- qrnn.predict(x = kern, parms = parms)

matplot(x, cbind(y, p), type=c("p", "l"), pch = 1, lwd = c(1, 2))
}
