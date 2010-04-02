\name{qrnn.cost}
\alias{qrnn.cost}
\title{
Huber norm approximation to the tilted absolute value cost function
}
\description{
Huber norm approximation to the tilted absolute value cost function
used to fit a QRNN model. Optional left censoring is supported.
}
\usage{
qrnn.cost(weights, x, y, n.hidden, tau, lower, eps, Th,
          Th.prime, penalty)
}
\arguments{
  \item{weights}{
      weight vector of length returned by \code{\link{qrnn.initialize}}.
  }
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
  }
  \item{y}{
      predictand column matrix with number of rows equal to the number of samples.
  }
  \item{n.hidden}{
      number of hidden nodes in the QRNN model.
  }
  \item{tau}{
      desired tau-quantile.
  }
  \item{lower}{
      left censoring point.
  }
  \item{eps}{
      epsilon value used in \code{\link{huber}} and related functions.
  }
  \item{Th}{
      hidden layer transfer function; use \code{\link{sigmoid}} for a nonlinear model and \code{\link{linear}} for a linear model.
  }
  \item{Th.prime}{
      derivative of the hidden layer transfer function \code{Th}.
  }
  \item{penalty}{
      weight penalty for weight decay regularization.
  }
}
\value{
    numeric value indicating tilted absolute value cost function, along with attribute containing vector with gradient information.
}
\seealso{
    \code{\link{qrnn.nlm}}, \code{\link{qrnn.fit}}
}
