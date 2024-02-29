\name{mcqrnn}
\alias{mcqrnn}
\alias{mcqrnn.fit}
\alias{mcqrnn.predict}
\title{
Monotone composite quantile regression neural network (MCQRNN) for simultaneous estimation of multiple non-crossing quantiles
}
\description{
High level wrapper functions for fitting and making predictions from a
monotone composite quantile regression neural network (MCQRNN) model for
multiple non-crossing regression quantiles (Cannon, 2018).

Uses \code{composite.stack} and monotonicity constraints in
\code{qrnn.fit} or \code{qrnn2.fit} to fit MCQRNN models with
one or two hidden layers. Note: \code{Th} must be a non-decreasing
function to guarantee non-crossing.

Following Tagasovska and Lopez-Paz (2019), it is also possible to estimate the
full quantile regression process by specifying a single integer value for 
\code{tau}. In this case, tau is the number of random samples used in the
stochastic estimation. It may be necessary to restart the optimization multiple
times from the previous weights and biases, in which case \code{init.range} can
be set to the weights values from the previously completed optimization run.
For large datasets, it is recommended that the \code{adam} method with an
appropriate \code{minibatch} size be used for optimization.

}
\usage{
mcqrnn.fit(x, y, n.hidden=2, n.hidden2=NULL, w=NULL,
           tau=c(0.1, 0.5, 0.9), iter.max=5000, n.trials=5,
           lower=-Inf, init.range=c(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
           monotone=NULL, eps.seq=2^seq(-8, -32, by=-4), Th=sigmoid,
           Th.prime=sigmoid.prime, penalty=0, n.errors.max=10,
           trace=TRUE, method=c("nlm", "adam"), scale.y=TRUE, ...)
mcqrnn.predict(x, parms, tau=NULL)
}
\arguments{
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
  }
  \item{y}{
      response column matrix with number of rows equal to the number of samples.
  }
  \item{n.hidden}{
      number of hidden nodes in the first hidden layer.
  }
   \item{n.hidden2}{
      number of hidden nodes in the second hidden layer; \code{NULL} fits a model with a single hidden layer.
  }
  \item{w}{
      if \code{tau} specifies a finite number of tau-quantiles, a vector of weights with length equal to the number of samples
      times the length of \code{tau}; see \code{composite.stack}. Otherwise, a vector of weights with length equal to the
      number of samples. \code{NULL} gives equal weight to each sample. 
  } 
  \item{tau}{
      desired tau-quantiles; \code{NULL} in \code{mcqrnn.predict} uses values from the original call to \code{mcqrnn.fit}.
      If \code{tau} is an integer, specifies the number of random samples used for stochastic estimation of the full quantile
      regression process.
  }
  \item{iter.max}{
      maximum number of iterations of the optimization algorithm.
  }
  \item{n.trials}{
      number of repeated trials used to avoid local minima.
  }
  \item{lower}{
      left censoring point.
  }
  \item{init.range}{
      initial weight range for input-hidden, hidden-hidden, and hidden-output weight matrices. If supplied with a list
      of weight matrices from a prior run of \code{mcqrnn.fit}, will restart model fitting with these values. 
  }
  \item{monotone}{
      column indices of covariates for which the monotonicity constraint should hold.
  }
  \item{eps.seq}{
      sequence of \code{eps} values for the finite smoothing algorithm.
  }
  \item{Th}{
      hidden layer transfer function; use \code{\link{sigmoid}}, \code{\link{elu}}, \code{\link{relu}},
      \code{\link{lrelu}}, \code{\link{softplus}}, or other non-decreasing function.
  }
  \item{Th.prime}{
      derivative of the hidden layer transfer function \code{Th}.
  }
  \item{penalty}{
      weight penalty for weight decay regularization.
  }
  \item{n.errors.max}{
      maximum number of \code{nlm} optimization failures allowed before quitting.
  }
  \item{trace}{
      logical variable indicating whether or not diagnostic messages are printed during optimization.
  }
  \item{method}{
      character string indicating which optimization algorithm to use when \code{n.hidden2 != NULL}.
  }
  \item{scale.y}{
      logical variable indicating whether \code{y} should be scaled to zero mean and unit standard deviation.
  }
  \item{\dots}{
      additional parameters passed to the \code{\link{nlm}} or \code{\link{adam}} optimization routines.
  }
  \item{parms}{
      list containing MCQRNN weight matrices and other parameters.
  }
}
\seealso{
    \code{\link{composite.stack}}, \code{\link{qrnn.fit}},
    \code{\link{qrnn2.fit}}, \code{\link{qrnn.predict}},
    \code{\link{qrnn2.predict}}, \code{\link{adam}}
}
\examples{
x <- as.matrix(iris[,"Petal.Length",drop=FALSE])
y <- as.matrix(iris[,"Petal.Width",drop=FALSE])

cases <- order(x)
x <- x[cases,,drop=FALSE]
y <- y[cases,,drop=FALSE]

set.seed(1)

## MCQRNN model w/ 2 hidden layers for simultaneous estimation of
## multiple non-crossing quantile functions
fit.mcqrnn <- mcqrnn.fit(x, y, tau=seq(0.1, 0.9, by=0.1),
                         n.hidden=2, n.hidden2=2, n.trials=1,
                         iter.max=500)
pred.mcqrnn <- mcqrnn.predict(x, fit.mcqrnn)

## Estimate the full quantile regression process by specifying
## the number of samples/random values of tau used in training

fit.full <- mcqrnn.fit(x, y, tau=1000L, n.hidden=3, n.hidden2=3,
                       n.trials=1, iter.max=300, eps.seq=1e-6,
                       method="adam", minibatch=64, print.level=100)
# Show how to initialize from previous weights
fit.full <- mcqrnn.fit(x, y, tau=1000L, n.hidden=3, n.hidden2=3,
                       n.trials=1, iter.max=300, eps.seq=1e-6,
                       method="adam", minibatch=64, print.level=100,
                       init.range=fit.full$weights)
pred.full <- mcqrnn.predict(x, fit.full, tau=seq(0.1, 0.9, by=0.1))

par(mfrow=c(1, 2))
matplot(x, pred.mcqrnn, col="blue", type="l")
points(x, y)
matplot(x, pred.full, col="blue", type="l")
points(x, y)

}
\references{
Cannon, A.J., 2011. Quantile regression neural networks: implementation
in R and application to precipitation downscaling. Computers & Geosciences,
37: 1277-1284. doi:10.1016/j.cageo.2010.07.005

Cannon, A.J., 2018. Non-crossing nonlinear regression quantiles by
monotone composite quantile regression neural network, with application
to rainfall extremes. Stochastic Environmental Research and Risk Assessment,
32(11): 3207-3225. doi:10.1007/s00477-018-1573-6

Tagasovska, N., D. Lopez-Paz, 2019. Single-model uncertainties for deep
learning. Advances in Neural Information Processing Systems, 32,
NeurIPS 2019. doi:10.48550/arXiv.1811.00908
}
