\name{qrnn.initialize}
\alias{qrnn.initialize}
\title{
Initialize a QRNN weight vector
}
\description{
Random initialization of the weight vector used during fitting of a QRNN model.
}
\usage{
qrnn.initialize(x, y, n.hidden, init.range=c(-0.5, 0.5, -0.5, 0.5))
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
  \item{init.range}{
      initial weight range for input-hidden and hidden-output weight matrices.
  }
}

