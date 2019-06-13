\name{qrnn-package}
\alias{qrnn-package}
\alias{qrnn}
\docType{package}
\title{Quantile Regression Neural Network}
\encoding{UTF-8}
\description{
This package implements the quantile regression neural network (QRNN)
(Taylor, 2000; Cannon, 2011; Cannon, 2018), which is a flexible nonlinear form
of quantile regression. The goal of quantile regression is to estimate
conditional quantiles of a response variable that depend on covariates
in some form of regression equation. The QRNN adopts the multi-layer
perceptron neural network architecture. The implementation follows from
previous work on the estimation of censored regression quantiles, thus
allowing predictions for mixed discrete-continuous variables like precipitation
(Friederichs and Hense, 2007). A differentiable approximation to the quantile
regression cost function is adopted so that a simplified form of the finite
smoothing algorithm (Chen, 2007) can be used to estimate model parameters.
This approximation can also be used to force the model to solve a standard
least squares regression problem or an expectile regression problem (Cannon, 2018).
Weight penalty regularization can be added to help avoid overfitting, and
ensemble models with bootstrap aggregation are also provided.

An optional monotone constraint can be invoked, which guarantees monotonic
non-decreasing behaviour of model outputs with respect to specified covariates
(Zhang, 1999). The input-hidden layer weight matrix can also be constrained
so that model relationships are strictly additive (see \code{\link{gam.style}};
Cannon, 2018). Borrowing strength by using a composite model for multiple
regression quantiles (Zou et al., 2008; Xu et al., 2017) is also possible
(see \code{\link{composite.stack}}). Weights can be applied to individual
cases (Jiang et al., 2012).

Applying the monotone constraint in combination with the composite model allows
one to simultaneously estimate multiple non-crossing quantiles (Cannon, 2018);
the resulting monotone composite QRNN (MCQRNN) is provided by the
\code{\link{mcqrnn.fit}} and \code{\link{mcqrnn.predict}} wrapper functions.
Examples for \code{\link{qrnn.fit}} and \code{\link{qrnn2.fit}} show how the
same functionality can be achieved using the low level \code{\link{composite.stack}}
and fitting functions.

QRNN models with a single layer of hidden nodes can be fitted using the
\code{\link{qrnn.fit}} function. Predictions from a fitted model are made using
the \code{\link{qrnn.predict}} function. Note: a single hidden layer
is usually sufficient for most modelling tasks. With added monotonicity
constraints, a second hidden layer may sometimes be beneficial
(Lang, 2005; Minin et al., 2010). QRNN models with two hidden layers are
available using the \code{\link{qrnn2.fit}} and
\code{\link{qrnn2.predict}} functions. For non-crossing quantiles, the
\code{\link{mcqrnn.fit}} wrapper also allows models with one or two hidden
layers to be fitted.

If models for multiple quantiles have been fitted, for example by
\code{\link{mcqrnn.fit}} or multiple calls to either \code{\link{qrnn.fit}}
or \code{\link{qrnn2.fit}}, the (experimental) \code{\link{dquantile}}
function and its companion functions are available to create proper
probability density, distribution, and quantile functions
(Quiñonero-Candela et al., 2006; Cannon, 2011). These can be useful for
assessing probabilistic calibration and evaluating model performance.

Finally, the function \code{\link{gam.style}} can be used to visualize and
investigate fitted covariate/response relationships from \code{qrnn.fit}
(Plate et al., 2000).

}
\details{
\tabular{ll}{
Package: \tab qrnn\cr
Type: \tab Package\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\references{
Cannon, A.J., 2011. Quantile regression neural networks: implementation
in R and application to precipitation downscaling. Computers & Geosciences,
37: 1277-1284. doi:10.1016/j.cageo.2010.07.005

Cannon, A.J., 2018. Non-crossing nonlinear regression quantiles by
monotone composite quantile regression neural network, with application
to rainfall extremes. Stochastic Environmental Research and Risk Assessment,
32(11): 3207-3225. doi:10.1007/s00477-018-1573-6

Chen, C., 2007. A finite smoothing algorithm for quantile regression.
Journal of Computational and Graphical Statistics, 16: 136-164.

Friederichs, P. and A. Hense, 2007. Statistical downscaling of extreme
precipitation events using censored quantile regression. Monthly Weather
Review, 135: 2365-2378. 

Jiang, X., J. Jiang, and X. Song, 2012. Oracle model selection for nonlinear
models based on weighted composite quantile regression. Statistica Sinica,
22(4): 1479-1506.

Lang, B., 2005. Monotonic multi-layer perceptron networks as universal
approximators. International Conference on Artificial Neural Networks,
Artificial Neural Networks: Formal Models and Their Applications-ICANN 2005,
pp. 31-37.

Minin, A., M. Velikova, B. Lang, and H. Daniels, 2010. Comparison of universal
approximators incorporating partial monotonicity by structure.
Neural Networks, 23(4): 471-475.

Plate, T., J. Bert, J. Grace, and P. Band, 2000. Visualizing the function
computed by a feedforward neural network. Neural Computation,
12(6): 1337-1354.

Quiñonero-Candela, J., C. Rasmussen, F. Sinz, O. Bousquet,
B. Scholkopf, 2006. Evaluating Predictive Uncertainty Challenge.
Lecture Notes in Artificial Intelligence, 3944: 1-27.

Taylor, J.W., 2000.  A quantile regression neural network approach to
estimating the conditional density of multiperiod returns. Journal of
Forecasting, 19(4): 299-311.

Xu, Q., K. Deng, C. Jiang, F. Sun, and X. Huang, 2017. Composite quantile
regression neural network with applications. Expert Systems with Applications,
76, 129-139.

Zhang, H. and Zhang, Z., 1999. Feedforward networks with monotone
constraints. In: International Joint Conference on Neural Networks,
vol. 3, p. 1820-1823. doi:10.1109/IJCNN.1999.832655

Zou, H. and M. Yuan, 2008. Composite quantile regression and the oracle model
selection theory. The Annals of Statistics, 1108-1126.
}
