\name{composite.stack}
\alias{composite.stack}
\title{
Reformat data matrices for composite quantile regression
}
\description{
Returns stacked \code{x} and \code{y} matrices and \code{tau} vector,
which can be passed to \code{qrnn.fit} to fit composite quantile
regression and composite QRNN models (Zou et al., 2008;
Xu et al., 2017). In combination with the partial monotonicity
constraints, stacking can be used to fit multiple non-crossing
quantile functions (see \code{\link{mcqrnn}}). More details
are provided in Cannon (2018).
}
\usage{
composite.stack(x, y, tau)
}
\arguments{
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to 
the number of variables.
  }
  \item{y}{
      response column matrix with number of rows equal to the number of samples.
  }
  \item{tau}{
      vector of tau-quantiles.
  }
}
\seealso{
    \code{\link{qrnn.fit}}, \code{\link{mcqrnn}}
}
\examples{
x <- as.matrix(iris[,"Petal.Length",drop=FALSE])
y <- as.matrix(iris[,"Petal.Width",drop=FALSE])

cases <- order(x)
x <- x[cases,,drop=FALSE]
y <- y[cases,,drop=FALSE]

tau <- seq(0.05, 0.95, by=0.05)
x.y.tau <- composite.stack(x, y, tau)
binary.tau <- dummy.code(as.factor(x.y.tau$tau))

set.seed(1)

# Composite QR
fit.cqr <- qrnn.fit(cbind(binary.tau, x.y.tau$x), x.y.tau$y,
                    tau=x.y.tau$tau, n.hidden=1, n.trials=1,
                    Th=linear, Th.prime=linear.prime)
pred.cqr <- matrix(qrnn.predict(cbind(binary.tau, x.y.tau$x), fit.cqr),
                   ncol=length(tau))
coef.cqr <- lm.fit(cbind(1, x), pred.cqr)$coef
colnames(coef.cqr) <- tau
print(coef.cqr)

# Composite QRNN
fit.cqrnn <- qrnn.fit(x.y.tau$x, x.y.tau$y, tau=x.y.tau$tau,
                      n.hidden=1, n.trials=1, Th=sigmoid,
                      Th.prime=sigmoid.prime)
pred.cqrnn <- qrnn.predict(x.y.tau$x, fit.cqrnn)
pred.cqrnn <- matrix(pred.cqrnn, ncol=length(tau), byrow=FALSE)

matplot(x, pred.cqrnn, col="red", type="l")
points(x, y, pch=20)

}

\references{
Cannon, A.J., 2018. Non-crossing nonlinear regression quantiles by monotone
composite quantile regression neural network, with application to rainfall
extremes. Stochastic Environmental Research and Risk Assessment
<https://dx.doi.org/10.1007/s00477-018-1573-6>. doi:10.1007/s00477-018-1573-6

Xu, Q., K. Deng, C. Jiang, F. Sun, and X. Huang, 2017. Composite quantile 
regression neural network with applications. Expert Systems with Applications,
76, 129-139.

Zou, H. and M. Yuan, 2008. Composite quantile regression and the oracle model
selection theory. The Annals of Statistics, 1108-1126.
}
