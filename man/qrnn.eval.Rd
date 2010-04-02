\name{qrnn.eval}
\alias{qrnn.eval}
\title{
Evaluate outputs from a QRNN model
}
\description{
Helper function used to evaluate quantiles from a fitted QRNN model.
}
\usage{
qrnn.eval(x, W1, W2, lower, eps, Th)
}
\arguments{
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
  }
  \item{W1}{
      input-hidden layer weight matrix.
  }
  \item{W2}{
      output-hidden layer weight matrix.
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
}
\value{
    a column matrix of predicted quantiles.
}
\seealso{
    \code{\link{qrnn.nlm}}, \code{\link{qrnn.predict}}
}
