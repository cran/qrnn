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
\code{tau}-quantile by minimizing a cost function based on smooth
Huber-norm approximations to the tilted absolute value and ramp functions.
Left censoring can be turned on by setting \code{lower} to a value
greater than \code{-Inf}. A simplified form of the finite smoothing
algorithm, in which the \code{\link{nlm}} optimization algorithm
is run with values of the \code{eps} approximation tolerance progressively
reduced in magnitude over the sequence \code{eps.seq}, is used to set the
QRNN weights and biases. Local minima of the cost function can be
avoided by setting \code{n.trials}, which controls the number of
repeated runs from different starting weights and biases, to a value
greater than one.

(Note: if \code{eps.seq} is set to a single, sufficiently large value and \code{tau}
is set to \code{0.5}, then the result will be a standard least squares
regression model. The same value of \code{eps.seq} and other values
of \code{tau} leads to expectile regression.)

The hidden layer transfer function \code{Th} and its derivative
\code{Th.prime} should be set to \code{\link{sigmoid}},
\code{\link{elu}}, or \code{\link{softplus}} and
\code{\link{sigmoid.prime}}, \code{\link{elu.prime}}, or \code{\link{softplus.prime}}
for a nonlinear model and to \code{\link{linear}} and
\code{\link{linear.prime}} for a linear model.

If invoked, the \code{monotone} argument enforces non-decreasing behaviour
between specified columns of \code{x} and model outputs. This holds if
\code{Th} and \code{To} are monotone non-decreasing functions. In this case,
the \code{exp} function is applied to the relevant weights following
initialization and during optimization; manual adjustment of
\code{init.weights} or \code{qrnn.initialize} may be needed due to
differences in scaling of the constrained and unconstrained weights.
Non-increasing behaviour can be forced by transforming the relevant
covariates, e.g., by reversing sign.

The \code{additive} argument sets relevant input-hidden layer weights
to zero, resulting in a purely additive model. Interactions between covariates
are thus suppressed, leading to a compromise in flexibility between
linear quantile regression and the quantile regression neural network.

Borrowing strength by using a composite model for multiple regression quantiles
is also possible (see \code{\link{composite.stack}}). Applying the monotone
constraint in combination with the composite model allows
one to simultaneously estimate multiple non-crossing quantiles;
the resulting monotone composite QRNN (MCQRNN) is demonstrated in
\code{\link{mcqrnn}}.

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

The \code{\link{gam.style}} function can be used to plot modified
generalized additive model effects plots, which are useful for visualizing
the modelled covariate-response relationships.

Note: values of \code{x} and \code{y} need not be standardized or
rescaled by the user. All variables are automatically scaled to zero
mean and unit standard deviation prior to fitting and parameters are
automatically rescaled by \code{\link{qrnn.predict}} and other prediction
functions. Values of \code{eps.seq} are relative to the residuals in
standard deviation units.
}
\usage{
qrnn.fit(x, y, n.hidden, w=NULL, tau=0.5, n.ensemble=1,
         iter.max=5000, n.trials=5, bag=FALSE, lower=-Inf,
         init.range=c(-0.5, 0.5, -0.5, 0.5), monotone=NULL,
         additive=FALSE, eps.seq=2^seq(-8, -32, by=-4),
         Th=sigmoid, Th.prime=sigmoid.prime, penalty=0,
         unpenalized=NULL, n.errors.max=10, trace=TRUE, ...)
}
\arguments{
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
      initial weight range for input-hidden and hidden-output weight matrices.
  }
  \item{monotone}{
      column indices of covariates for which the monotonicity constraint should hold.
  }
  \item{additive}{
      force additive relationships.
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
  \item{\dots}{
      additional parameters passed to the \code{\link{nlm}} optimization routine.
  }
}
\value{
    a list containing elements
    \item{weights}{a list containing fitted weight matrices}
    \item{lower}{left censoring point}
    \item{eps.seq}{sequence of \code{eps} values for the finite smoothing algorithm}
    \item{tau}{desired tau-quantile(s)}
    \item{Th}{hidden layer transfer function}
    \item{x.center}{vector of column means for \code{x}}
    \item{x.scale}{vector of column standard deviations for \code{x}}
    \item{y.center}{vector of column means for \code{y}}
    \item{y.scale}{vector of column standard deviations for \code{y}}
    \item{monotone}{column indices indicating covariate monotonicity constraints.}
    \item{additive}{force additive relationships.}
}
\seealso{
    \code{\link{qrnn.predict}}, \code{\link{qrnn.cost}}, \code{\link{composite.stack}}, \code{\link{mcqrnn}}, \code{\link{gam.style}}
}
\examples{
x <- as.matrix(iris[,"Petal.Length",drop=FALSE])
y <- as.matrix(iris[,"Petal.Width",drop=FALSE])

cases <- order(x)
x <- x[cases,,drop=FALSE]
y <- y[cases,,drop=FALSE]

tau <- c(0.05, 0.5, 0.95)
 
set.seed(1)

## QRNN models for conditional 5th, 50th, and 95th percentiles
w <- p <- vector("list", length(tau))
for(i in seq_along(tau)){
    w[[i]] <- qrnn.fit(x=x, y=y, n.hidden=3, tau=tau[i],
                       iter.max=200, n.trials=1)
    p[[i]] <- qrnn.predict(x, w[[i]])
}

## Monotone composite QRNN (MCQRNN) for simultaneous estimation of
## multiple non-crossing quantile functions
x.y.tau <- composite.stack(x, y, tau)
fit.mcqrnn <- qrnn.fit(cbind(x.y.tau$tau, x.y.tau$x), x.y.tau$y,
                       tau=x.y.tau$tau, n.hidden=3, n.trials=1,
                       iter.max=500, monotone=1)
pred.mcqrnn <- matrix(qrnn.predict(cbind(x.y.tau$tau, x.y.tau$x),
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

Cannon, A.J., 2018. Non-crossing nonlinear regression quantiles by
monotone composite quantile regression neural network, with application
to rainfall extremes. Stochastic Environmental Research and Risk Assessment,
32(11): 3207-3225. doi:10.1007/s00477-018-1573-6
}
