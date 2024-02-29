### Description

The [qrnn](https://cran.r-project.org/package=qrnn) package for R implements
the quantile regression neural network (QRNN) (Taylor, 2000; Cannon, 2011;
Cannon, 2018), which is a flexible nonlinear form of quantile regression.
While low level modelling functions are available, it is recommended that the
`mcqrnn.fit` and `mcqrnn.predict` wrappers be used for most applications. More
information is provided below.

The goal of quantile regression is to estimate conditional quantiles of a
response variable that depend on covariates in some form of regression
equation. The QRNN adopts the multi-layer perceptron neural network
architecture. The implementation follows from previous work on the estimation
of censored regression quantiles, thus allowing predictions for mixed
discrete-continuous variables like precipitation (Friederichs and Hense, 2007).
A differentiable approximation to the quantile regression cost function is
adopted so that a simplified form of the finite smoothing algorithm
(Chen, 2007) can be used to estimate model parameters. This approximation can
also be used to force the model to solve a standard least squares regression
problem or an expectile regression problem (Cannon, 2018). Weight penalty
regularization can be added to help avoid overfitting, and ensemble models
with bootstrap aggregation are also provided.

An optional monotone constraint can be invoked, which guarantees monotonic
non-decreasing behaviour of model outputs with respect to specified covariates
(Zhang, 1999). The input-hidden layer weight matrix can also be constrained
so that model relationships are strictly additive (see `gam.style`;
Cannon, 2018). Borrowing strength by using a composite model for multiple
regression quantiles (Zou et al., 2008; Xu et al., 2017) is also possible
(see `composite.stack`). Weights can be applied to individual cases
(Jiang et al., 2012).

Applying the monotone constraint in combination with the composite model allows
one to simultaneously estimate multiple non-crossing quantiles (Cannon, 2018);
the resulting monotone composite QRNN (MCQRNN) is provided by the
`mcqrnn.fit` and `mcqrnn.predict` wrapper functions. Examples for `qrnn.fit`
and `qrnn2.fit` show how the same functionality can be achieved using the low
level `composite.stack` and fitting functions.

QRNN models with a single layer of hidden nodes can be fitted using the
`qrnn.fit` function. Predictions from a fitted model are made using
the `qrnn.predict` function. The function `gam.style` can be used to visualize
and investigate fitted covariate/response relationships from `qrnn.fit`
(Plate et al., 2000). Note: a single hidden layer is usually sufficient
for most modelling tasks. With added monotonicity constraints, a second hidden
layer may sometimes be beneficial (Lang, 2005; Minin et al., 2010). QRNN models
with two hidden layers are available using the `qrnn2.fit` and `qrnn2.predict`
functions. For non-crossing quantiles, the `mcqrnn.fit` and `mcqrnn.predict`
wrappers also allow models with one or two hidden layers to be fitted and
predictions to be made from the fitted models.

In general, `mcqrnn.fit` offers a convenient, single function for fitting
multiple quantiles simultaneously. Note, however, that default settings in
mcqrnn.fit and other model fitting functions are not optimized for general
speed, memory efficiency, or accuracy and should be adjusted for a particular
regression problem as needed. In particular, the approximation to the quantile
regression cost function `eps.seq`, the number of trials `n.trials`, and number
of iterations `iter.max` can all influence fitting speed (and accuracy), as can
changing the optimization algorithm via `method`. Non-crossing quantiles are 
implemented by stacking multiple copies of the `x` and `y` data, one copy per
value of `tau`. Depending on the dataset size, this can lead to large matrices
being passed to the optimization routine. In the `adam` adaptive stochastic
gradient descent method, the `minibatch` size can be adjusted to help offset
this cost. Model complexity is determined via the number of hidden nodes, 
`n.hidden` and `n.hidden2`, as well as the optional weight penalty `penalty`;
values of these hyperparameters are crucial to obtaining a well performing
model.

When using `mcqrnn.fit`, it is also possible to estimate the full quantile
regression process by specifying a single integer value for `tau`. In this
case, `tau` is the number of random samples used in the stochastic estimation.
For more information, see Tagasovska and Lopez-Paz (2019). It may be necessary
to restart the optimization multiple times from the previous weights and biases,
in which case `init.range` can be set to the `weights` values from the
previously completed optimization run. For large datasets, it is recommended
that the `adam` method with an integer `tau` and an appropriate `minibatch`
size be used for optimization.

If models for multiple quantiles have been fitted, for example by
`mcqrnn.fit` or multiple calls to either `qrnn.fit` or `qrnn2.fit`, the
(experimental) `dquantile` function and its companion functions are available to
create proper probability density, distribution, and quantile functions
(Quiñonero-Candela et al., 2006; Cannon, 2011). Alternative distribution,
quantile, and random variate functions based on the Nadaraya-Watson estimator
(Passow and Donner, 2020) are also available in `[p,q,r]quantile.nw`. These can
be useful for assessing probabilistic calibration and evaluating model
performance.

Note: the user cannot easily change the output layer transfer function
to be different than `hramp`, which provides either the identity function or a
ramp function to accommodate optional left censoring. Some applications, for
example fitting smoothed binary quantile regression models for a binary target
variable (Kordas, 2006), require an alternative like the logistic sigmoid.
While not straightforward, it is possible to change the output layer transfer
function by switching off `scale.y` in the call to the fitting
function and reassigning `hramp` and `hramp.prime` as follows:

```
library(qrnn)

# Use the logistic sigmoid as the output layer transfer function
To.logistic <- function(x, lower, eps) 0.5 + 0.5*tanh(x/2)
environment(To.logistic) <- asNamespace("qrnn")
assignInNamespace("hramp", To.logistic, ns="qrnn")

# Change the derivative of the output layer transfer function
To.logistic.prime <- function(x, lower, eps) 0.25/(cosh(x/2)^2)
environment(To.logistic.prime) <- asNamespace("qrnn")
assignInNamespace("hramp.prime", To.logistic.prime, ns="qrnn")

```

### References

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

Kordas, G., 2006. Smoothed binary regression quantiles. Journal of Applied
Econometrics, 21(3): 387-407.

Lang, B., 2005. Monotonic multi-layer perceptron networks as universal
approximators. International Conference on Artificial Neural Networks,
Artificial Neural Networks: Formal Models and Their Applications-ICANN 2005,
pp. 31-37.

Minin, A., M. Velikova, B. Lang, and H. Daniels, 2010. Comparison of universal
approximators incorporating partial monotonicity by structure.
Neural Networks, 23(4): 471-475.

Passow, C., R.V. Donner, 2020. Regression-based distribution mapping for
bias correction of climate model outputs using linear quantile regression.
Stochastic Environmental Research and Risk Assessment, 34: 87-102.

Plate, T., J. Bert, J. Grace, and P. Band, 2000. Visualizing the function
computed by a feedforward neural network. Neural Computation,
12(6): 1337-1354.

Quiñonero-Candela, J., C. Rasmussen, F. Sinz, O. Bousquet,
B. Scholkopf, 2006. Evaluating Predictive Uncertainty Challenge.
Lecture Notes in Artificial Intelligence, 3944: 1-27.

Tagasovska, N., D. Lopez-Paz, 2019. Single-model uncertainties for deep
learning. Advances in Neural Information Processing Systems, 32,
NeurIPS 2019. doi:10.48550/arXiv.1811.00908

Taylor, J.W., 2000. A quantile regression neural network approach to
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
