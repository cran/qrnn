\name{qrnn2}
\alias{qrnn2.fit}
\alias{qrnn2.predict}
\title{
Fit and make predictions from QRNN models with two hidden layers
}
\description{
Functions used to fit and make predictions from QRNN models with two hidden layers.
Note: \code{Th} must be a non-decreasing function if \code{monotone != NULL}.
}
\usage{
qrnn2.fit(x, y, n.hidden=2, n.hidden2=2, w=NULL, tau=0.5,
          n.ensemble=1, iter.max=5000, n.trials=5, bag=FALSE,
          lower=-Inf, init.range=c(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
          monotone=NULL, eps.seq=2^(-8:-32), Th=sigmoid,
          Th.prime=sigmoid.prime, penalty=0, unpenalized=NULL,
          n.errors.max=10, trace=TRUE, method=c("nlm", "adam"),
          ...)
qrnn2.predict(x, parms)
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
      number of hidden nodes in the second hidden layer.
  }
  \item{w}{
      vector of weights with length equal to the number of samples;
      \code{NULL} gives equal weight to each sample.
  } 
  \item{tau}{
      desired tau-quantile(s).
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
  \item{init.range}{
      initial weight range for input-hidden, hidden-hidden, and hidden-output weight matrices.
  }
  \item{monotone}{
      column indices of covariates for which the monotonicity constraint should hold.
  }
  \item{eps.seq}{
      sequence of \code{eps} values for the finite smoothing algorithm.
  }
  \item{Th}{
      hidden layer transfer function; use \code{\link{sigmoid}}, \code{\link{elu}}, or \code{\link{softplus}} for a nonlinear model and \code{\link{linear}} for a linear model.
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
  \item{n.errors.max}{
      maximum number of \code{nlm} optimization failures allowed before quitting.
  }
  \item{trace}{
      logical variable indicating whether or not diagnostic messages are printed during optimization.
  }
  \item{method}{
      character string indicating which optimization algorithm to use.
  }
  \item{\dots}{
      additional parameters passed to the \code{\link{nlm}} or \code{\link{adam}} optimization routines.
  }
  \item{parms}{
      list containing QRNN weight matrices and other parameters from \code{\link{qrnn2.fit}}.
  }
}
\seealso{
    \code{\link{qrnn.fit}}, \code{\link{qrnn.predict}},
    \code{\link{qrnn.cost}}, \code{\link{composite.stack}},
    \code{\link{mcqrnn}}, \code{\link{adam}}
}
\examples{
x <- as.matrix(iris[,"Petal.Length",drop=FALSE])
y <- as.matrix(iris[,"Petal.Width",drop=FALSE])

cases <- order(x)
x <- x[cases,,drop=FALSE]
y <- y[cases,,drop=FALSE]

tau <- c(0.05, 0.5, 0.95)
 
set.seed(1)

## QRNN models w/ 2 hidden layers (tau=0.05, 0.50, 0.95)
w <- p <- vector("list", length(tau))
for(i in seq_along(tau)){
    w[[i]] <- qrnn2.fit(x=x, y=y, n.hidden=3, n.hidden2=3,
                       tau=tau[i], iter.max=200, n.trials=1)
    p[[i]] <- qrnn2.predict(x, w[[i]])
}

## MCQRNN model w/ 2 hidden layers for simultaneous estimation of
## multiple non-crossing quantile functions
x.y.tau <- composite.stack(x, y, tau)
fit.mcqrnn <- qrnn2.fit(cbind(x.y.tau$tau, x.y.tau$x), x.y.tau$y,
                        tau=x.y.tau$tau, n.hidden=3, n.hidden2=3,
                        n.trials=1, iter.max=500, monotone=1)
pred.mcqrnn <- matrix(qrnn2.predict(cbind(x.y.tau$tau, x.y.tau$x),
                      fit.mcqrnn), ncol=length(tau))

par(mfrow=c(1, 2))
matplot(x, matrix(unlist(p), nrow=nrow(x), ncol=length(p)), col="red",
        type="l")
points(x, y)
matplot(x, pred.mcqrnn, col="blue", type="l")
points(x, y)

}
\references{
Cannon, A.J., 2011. Quantile regression neural networks: implementation
in R and application to precipitation downscaling. Computers & Geosciences,
37: 1277-1284. doi:10.1016/j.cageo.2010.07.005

Cannon, A.J., 2018. Non-crossing nonlinear regression quantiles by monotone
composite quantile regression neural network, with application to rainfall
extremes. Stochastic Environmental Research and Risk Assessment
<https://dx.doi.org/10.1007/s00477-018-1573-6>. doi:10.1007/s00477-018-1573-6
}
