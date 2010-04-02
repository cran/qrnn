\name{qrnn.predict}
\alias{qrnn.predict}
\title{
Evaluate quantiles from trained QRNN model
}
\description{
Evaluate a fitted QRNN model or ensemble of models, resulting in a list
containing the predicted quantiles.
}
\usage{
qrnn.predict(x, parms)
}
\arguments{
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
  }
  \item{parms}{
      list containing QRNN input-hidden and hidden-output layer weight matrices and other parameters from \code{\link{qrnn.fit}}.
  }
}
\value{
    a list with number of elements equal to that of \code{parms}, each containing a column matrix of predicted quantiles.
}
\seealso{
    \code{\link{qrnn.fit}}
}
\examples{
data(sinc)
x <- sinc$x
y <- sinc$y
y[y < -0.1] <- -0.1

parms <- qrnn.fit(x = x, y = y, n.hidden = 5, tau = 0.5, lower = -0.1,
                  iter.max = 1000, n.trials = 2)
p <- qrnn.predict(x = x, parms = parms)

matplot(x, cbind(y, p), type=c("p", "l"), pch = 1, lwd = c(1, 2))
}
