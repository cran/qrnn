\name{qrnn.cost}
\alias{qrnn.cost}
\title{
Smooth approximation to the tilted absolute value cost function
}
\description{
Smooth approximation to the tilted absolute value cost function
used to fit a QRNN model. Optional left censoring, monotone constraints,
and additive constraints are supported.
}
\usage{
qrnn.cost(weights, x, y, n.hidden, w, tau, lower, monotone,
          additive, eps, Th, Th.prime, penalty, unpenalized)
}
\arguments{
  \item{weights}{
      weight vector of length returned by \code{\link{qrnn.initialize}}.
  }
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
  }
  \item{y}{
      response column matrix with number of rows equal to the number of samples.
  }
  \item{n.hidden}{
      number of hidden nodes in the QRNN model.
  }
  \item{w}{
      vector of weights with length equal to the number of samples;
      \code{NULL} gives equal weight to each sample.
  }
  \item{tau}{
      desired tau-quantile.
  }
  \item{lower}{
      left censoring point.
  }
  \item{monotone}{
      column indices of covariates for which the monotonicity constraint should hold.
  }
  \item{additive}{
      force additive relationships.
  }
  \item{eps}{
      epsilon value used in the approximation functions.
  }
  \item{Th}{
      hidden layer transfer function; use \code{\link{sigmoid}}, \code{\link{elu}}, \code{\link{relu}}, \code{\link{lrelu}}, 
      \code{\link{softplus}}, or other non-decreasing function for a nonlinear model and \code{\link{linear}} for a linear model.
  }
  \item{Th.prime}{
      derivative of the hidden layer transfer function \code{Th}.
  }
  \item{penalty}{
      weight penalty for weight decay regularization.
  }
  \item{unpenalized}{
      column indices of covariates for which the weight penalty should not be applied to input-hidden layer weights.
  }
}
\value{
    numeric value indicating tilted absolute value cost function, along with attribute containing vector with gradient information.
}
\seealso{
    \code{\link{qrnn.fit}}
}

