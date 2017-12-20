composite.stack <- function(x, y, tau){
    n <- nrow(x)
    x <- kronecker(rep(1, length(tau)), x)
    y <- kronecker(rep(1, length(tau)), y)
    tau <- rep(tau, each=n)
    list(x=x, y=y, tau=tau)
}

