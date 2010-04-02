\name{qrnn.initialize}
\alias{qrnn.initialize}
\title{
Initialize a QRNN weight vector
}
\description{
Random initialization of the weight vector used during fitting of a QRNN model.
}
\usage{
qrnn.initialize(x, y, n.hidden)
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
}
\seealso{
    \code{\link{qrnn.reshape}}
}
