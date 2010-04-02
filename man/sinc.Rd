\name{sinc}
\alias{sinc}
\docType{data}
\title{
Noisy sinc dataset
}
\description{
Noisy sinc dataset generated via
\code{rnorm(200, mean=sin(x)/x, sd=0.2*abs(sin(x)/x))} with
\code{x <- seq(-3*pi, 3*pi, length=200)}. Includes the theoretical
5th, 50th, and 95th percentiles.
}
\examples{
data(sinc)
plot(sinc$x, sinc$y, ylim = range(pretty(c(sinc$y, sinc$tau))))
matlines(sinc$x, sinc$tau, lwd = 2)
}
\keyword{datasets}
