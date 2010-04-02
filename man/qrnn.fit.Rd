\name{qrnn.fit}
\alias{qrnn.fit}
\title{
Main function used to fit a QRNN model or ensemble of QRNN models
}
\description{
Function used to fit a QRNN model or ensemble of QRNN models.
}
\details{
Fit a censored quantile regression neural network model for the
\code{tau}-quantile by minimizing a cost function based on the Huber
norm approximation to the tilted absolute value and ramp functions.
Left censoring can be turned on by setting \code{lower} to a value
greater than \code{-Inf}. A simplified form of the finite smoothing
algorithm, in which the \code{\link{nlm}} optimization algorithm
is run with values of the Huber norm \code{eps} parameter progressively
reduced in magnitude over the sequence \code{eps.seq}, is used to set the
QRNN weights and biases. Local minima of the cost function can be
avoided by setting \code{n.trials}, which controls the number of
repeated runs from different starting weights and biases, to a value
greater than one.

The hidden layer transfer function \code{Th} and its derivative
\code{Th.prime} should be set to \code{\link{sigmoid}} and
\code{\link{sigmoid.prime}} for a nonlinear model and to
\code{\link{linear}} and \code{\link{linear.prime}} for a linear
model.

In the linear case, model complexity does not depend on the number
of hidden nodes; the value of \code{n.hidden} is ignored and is instead
set to one internally. In the nonlinear case, \code{n.hidden}
controls the overall complexity of the model. As an added means of
avoiding overfitting, weight penalty regularization for the magnitude
of the input-hidden layer weights (excluding biases) can be applied
by setting \code{penalty} to a nonzero value. (For the linear model,
this penalizes both input-hidden and hidden-output layer weights,
leading to a quantile ridge regression model. In this case, kernel
quantile ridge regression can be performed with the aid of the
\code{\link{qrnn.rbf}} function.) Finally, if the \code{bag} argument
is set to \code{TRUE}, models are trained on bootstrapped \code{x} and
\code{y} sample pairs; bootstrap aggregation (bagging) can be turned
on by setting \code{n.ensemble} to a value greater than one. Averaging
over an ensemble of bagged models will also tend to alleviate
overfitting.

Note: values of \code{x} and \code{y} need not be standardized or
rescaled by the user. All variables are automatically scaled to zero
mean and unit standard deviation prior to fitting and parameters are
automatically rescaled by \code{\link{qrnn.predict}}. Values of
\code{eps.seq} are relative to the residuals in standard deviation
units.
}
\usage{
qrnn.fit(x, y, n.hidden, tau = 0.5, n.ensemble = 1,
         iter.max = 5000, n.trials = 5, bag = FALSE,
         lower = -Inf, eps.seq = 2^(-8:-32), Th = sigmoid,
         Th.prime = sigmoid.prime, penalty = 0, ...)
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
  \item{n.ensemble}{
     number of ensemble members to fit.
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
    \item{weights}{a list containing fitted weight matrices}
    \item{lower}{left censoring point}
    \item{eps.seq}{sequence of \code{eps} values for the finite smoothing algorithm}
    \item{tau}{desired tau-quantile}
    \item{Th}{hidden layer transfer function}
    \item{x.center}{vector of column means for \code{x}}
    \item{x.scale}{vector of column standard deviations for \code{x}}
    \item{y.center}{vector of column means for \code{y}}
    \item{y.scale}{vector of column standard deviations for \code{y}}
}
\seealso{
    \code{\link{qrnn.predict}}, \code{\link{qrnn.nlm}}, \code{\link{qrnn.cost}}
}
\examples{
data(sinc)
x <- sinc$x
y <- sinc$y
q <- sinc$tau
probs <- c(0.05, 0.50, 0.95)
 
## Fit QRNN models for 5th, 50th, and 95th percentiles

w <- p <- list()
for(i in seq_along(probs)){
    w[[i]] <- qrnn.fit(x = x, y = y, n.hidden = 5, tau = probs[i],
                       iter.max = 1000, n.trials = 2)
    p[[i]] <- qrnn.predict(x, w[[i]])
}

plot(x, y, ylim = range(pretty(c(y, q))))
matlines(x, q, lwd = 2)
matlines(x, matrix(unlist(p), nrow = nrow(x), ncol = length(p)))

}
