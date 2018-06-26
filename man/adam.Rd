\name{adam}
\alias{adam}
\title{
Adaptive stochastic gradient descent optimization algorithm (Adam)
}
\description{
From Kingma and Ba (2015): "We introduce Adam, an algorithm for first-order 
gradient-based optimization of stochastic objective functions, based on 
adaptive estimates of lower-order moments. The method is straightforward to 
implement, is computationally efficient, has little memory requirements, is 
invariant to diagonal rescaling of the gradients, and is well suited for 
problems that are large in terms of data and/or parameters. The method is 
also appropriate for non-stationary objectives and problems with very noisy 
and/or sparse gradients. The hyper-parameters have intuitive interpretations 
and typically require little tuning. Some connections to related algorithms, 
on which Adam was inspired, are discussed. We also analyze the theoretical 
convergence properties of the algorithm and provide a regret bound on the 
convergence rate that is comparable to the best known results under the 
online convex optimization framework. Empirical results demonstrate that Adam 
works well in practice and compares favorably to other stochastic 
optimization methods. Finally, we discuss AdaMax, a variant of Adam based on 
the infinity norm."
}
\usage{
adam(f, p, x, y, w, tau, ..., iterlim=5000, iterbreak=iterlim,
     alpha=0.01, minibatch=nrow(x), beta1=0.9, beta2=0.999,
     epsilon=1e-8, print.level=10)

}
\arguments{
  \item{f}{
the function to be minimized, including gradient information
contained in the \code{gradient} attribute.
}
  \item{p}{
the starting parameters for the minimization.
}
  \item{x}{
covariate matrix with number of rows equal to the number of
samples and number of columns equal to the number of variables.
}
  \item{y}{
response column matrix with number of rows equal to the number
of samples.
}
  \item{w}{
vector of weights with length equal to the number of samples.
}
  \item{tau}{
vector of desired tau-quantile(s) with length equal to the
number of samples.
}
  \item{\dots}{
      additional parameters passed to the \code{f} cost function.
}
  \item{iterlim}{
the maximum number of iterations before the optimization is stopped.
}
  \item{iterbreak}{
the maximum number of iterations without progress before the
optimization is stopped.
}
  \item{alpha}{
size of the learning rate.
}
  \item{minibatch}{
number of samples in each minibatch.
}
  \item{beta1}{
controls the exponential decay rate used to scale the biased first
moment estimate.
}
  \item{beta2}{
controls the exponential decay rate used to scale the biased second
raw moment estimate.
}
  \item{epsilon}{
smoothing term to avoid division by zero.
}
  \item{print.level}{
the level of printing which is done during optimization. A value of
\code{0} suppresses any progress reporting, whereas positive values
report the value of \code{f} every \code{print.level} iterations.
}
}
\value{
A list with elements:
  \item{estimate}{The best set of parameters found.}
  \item{minimum}{The value of \code{f} corresponding to \code{estimate}.}
}
\references{
Kingma, D.P. and J. Ba, 2015. Adam: A method for stochastic optimization.
The International Conference on Learning Representations (ICLR) 2015.
http://arxiv.org/abs/1412.6980
}
