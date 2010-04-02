\name{qrnn.nlm}
\alias{qrnn.nlm}
\title{
Fit a QRNN model via nlm()
}
\description{
Helper function used to fit a QRNN model via the nlm() function
and a variant of the finite smoothing algorithm.
}
\usage{
qrnn.nlm(x, y, n.hidden, tau, iter.max, n.trials, bag,
         lower, eps.seq, Th, Th.prime, penalty, ...)
}
\arguments{
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
  \item{iter.max}{
      maximum number of iterations of the optimization algorithm.
  }
  \item{n.trials}{
      number of repeated trials used to avoid local minima.
  }
  \item{bag}{
      logical variable indicating whether or not bootstrap aggregation (bagging) should be used.
  }
  \item{lower}{
      left censoring point.
  }
  \item{eps.seq}{
      sequence of \code{eps} values for the finite smoothing algorithm; used in \code{\link{huber}} and related functions.
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
  \item{\dots}{
      additional parameters passed to the \code{\link{nlm}} optimization routine.
  }
}
\value{
    a list containing elements
    \item{W1}{matrix of optimized input-hidden layer weights.}
    \item{W2}{matrix of optimized hidden-output layer weights.}
}
\seealso{
    \code{\link{qrnn.cost}}, \code{\link{qrnn.fit}}, \code{\link{qrnn.eval}}
}
