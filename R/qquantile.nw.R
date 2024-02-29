###############################################################################
# Quantile, cumulative distribution, and random variate functions based on
# interpolation of [tau, quant] pairs using the Nadaraya-Watson estimator
# with a beta kernel (bandwith = h):
# Passow, C., R.V. Donner, 2020. Regression-based distribution mapping for
# bias correction of climate model outputs using linear quantile regression.
#  Stochastic Environmental Research and Risk Assessment, 34:87-102.
#  doi:10.1007/s00477-019-01750-7

qquantile.nw <- function(p, tau, quant, h=0.001){
# Quantile function based on [tau, quant] pairs
    K <- function(p, tau, h){
        ((p^(tau/h))*(1-p)^((1-tau)/h))/
        beta(tau/h + 1, (1-tau)/h + 1)
    }
    q <- sum(K(p, tau, h)*quant)/sum(K(p, tau, h))
    q
}

pquantile.nw <- function(q, tau, quant, h=0.001, ...){
# Cumulative distribution function based on [tau, quant] pairs
    func <- function(p, q, tau, quant, h){
        qq <- qquantile.nw(p, tau=tau, quant=quant, h=h)
           q-qq
    }
    p <- uniroot(f=func, q=q, lower=min(tau), upper=max(tau), tau=tau,
                 quant=quant, h=h, ...)$root
    p
}

rquantile.nw <- function(n, tau, quant, h=0.001){
# Random variate function based on [tau, quant] pairs
    sapply(runif(n), qquantile.nw, tau=tau, quant=quant, h=h)
}

################################################################################