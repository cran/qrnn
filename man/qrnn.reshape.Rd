\name{qrnn.reshape}
\alias{qrnn.reshape}
\title{
Reshape a QRNN weight vector
}
\description{
Reshapes a weight vector used during fitting of a QRNN model into input-hidden and hidden-output layer weight matrices.
}
\usage{
qrnn.reshape(x, y, weights, n.hidden)
}
\arguments{
  \item{x}{
      covariate matrix with number of rows equal to the number of samples and number of columns equal to the number of variables.
  }
  \item{y}{
      predictand column matrix with number of rows equal to the number of samples.
  }
  \item{weights}{
      weight vector of length returned by \code{\link{qrnn.initialize}}.
  }
  \item{n.hidden}{
      number of hidden nodes in the QRNN model.
  }
}
\seealso{
    \code{\link{qrnn.initialize}}
}
