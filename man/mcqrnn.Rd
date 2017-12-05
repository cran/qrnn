\name{mcqrnn}
\alias{mcqrnn}
\title{
Monotone composite quantile regression neural network (MCQRNN)
}
\description{
Use \code{composite.stack} and monotonicity constraints in \code{qrnn.fit}
to fit a monotone composite quantile regression model (MCQRNN) for
multiple non-crossing quantiles (Cannon, 2017).
}
\seealso{
    \code{\link{composite.stack}}, \code{\link{qrnn.fit}}
}
\examples{
x <- as.matrix(iris[,"Petal.Length",drop=FALSE])
y <- as.matrix(iris[,"Petal.Width",drop=FALSE])

cases <- order(x)
x <- x[cases,,drop=FALSE]
y <- y[cases,,drop=FALSE]

tau <- seq(0.05, 0.95, by=0.05)
x.y.tau <- composite.stack(x, y, tau)

set.seed(1)

## Monotone composite QRNN (MCQRNN) for simultaneous estimation of
## multiple non-crossing quantile functions
fit.mcqrnn <- qrnn.fit(cbind(x.y.tau$tau, x.y.tau$x), x.y.tau$y,
                       tau=x.y.tau$tau, n.hidden=3, n.trials=1,
                       iter.max=400, monotone=1)

## Add partial monotonicity constraint on the x covariate variable
fit.mcqrnn.mx <-  qrnn.fit(cbind(x.y.tau$tau, x.y.tau$x), x.y.tau$y,
                           tau=x.y.tau$tau, n.hidden=3, n.trials=1,
                           iter.max=400, monotone=c(1, 2))

pred.mcqrnn <- matrix(qrnn.predict(cbind(x.y.tau$tau, x.y.tau$x),
                      fit.mcqrnn), ncol=length(tau))
pred.mcqrnn.mx <- matrix(qrnn.predict(cbind(x.y.tau$tau, x.y.tau$x),
                         fit.mcqrnn.mx), ncol=length(tau))

par(mfrow=c(1, 2))
matplot(x, pred.mcqrnn, col="red", type="l")
points(x, y, pch=20)
matplot(x, pred.mcqrnn.mx, col="blue", type="l")
points(x, y, pch=20)
}

\references{
Cannon, A.J., 2017. Non-crossing nonlinear regression quantiles by monotone
composite quantile regression neural network, with application to rainfall
extremes. EarthArXiv <https://eartharxiv.org/wg7sn>. doi:10.17605/OSF.IO/WG7SN
}
