qrnn.fit <-
function(x, y, n.hidden, w=NULL, tau=0.5, n.ensemble=1, iter.max=5000,
         n.trials=5, bag=FALSE, lower=-Inf, init.range=c(-0.5, 0.5, -0.5, 0.5),
         monotone=NULL, additive=FALSE, eps.seq=2^(-8:-32), Th=sigmoid,
         Th.prime=sigmoid.prime, penalty=0, unpenalized=NULL,
         n.errors.max=10, trace=TRUE, ...)

{
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    if (!is.matrix(y)) stop("\"y\" must be a matrix")
    if (ncol(y) != 1) stop("\"y\" must be univariate")
    if (any((tau > 1) | (tau < 0))) stop("invalid \"tau\"")
    if (!identical(Th, linear) && missing(n.hidden)) stop("must specify \"n.hidden\"")
    if (identical(Th, linear)) n.hidden <- 1
    is.whole <- function(x, tol = .Machine$double.eps^0.5) 
        abs(x - round(x)) < tol
    if (additive && !is.whole(n.hidden/ncol(x))) stop("\"n.hidden\" must be an integer multiple of \"ncol(x)\" when \"additive=TRUE\"")
    if(is.null(w)) w <- rep(1/nrow(y), nrow(y))
    if (any(w < 0)) stop("invalid \"w\"")
    x <- scale(x)
    x.center <- attr(x, "scaled:center")
    x.scale <- attr(x, "scaled:scale")
    y <- scale(y)
    y.center <- attr(y, "scaled:center")
    y.scale <- attr(y, "scaled:scale")
    lower.scaled <- (lower-y.center)/y.scale
    if(additive)
        additive <- gam.mask(x, n.hidden)
    weights <- list()
    if(trace) cat("tau =", unique(tau), "\n", sep=" ")
    for (i in seq(n.ensemble)){
        if(trace) cat(i, "/", n.ensemble, "\n", sep="")
        w.tmp <- NA
        class(w.tmp) <- "try-error"
        n.errors <- 0
        while (inherits(w.tmp, "try-error")) {
            w.tmp <- try(qrnn.nlm(x, y, n.hidden, w, tau, iter.max,
                                  n.trials, bag, lower.scaled,
                                  init.range, monotone, additive, eps.seq,
                                  Th, Th.prime, penalty, unpenalized,
                                  trace, ...),
                        silent = TRUE)
            n.errors <- n.errors + 1
            if (n.errors > n.errors.max) stop("nlm optimization failed")
        }
        weights[[i]] <- w.tmp
    }
    if(trace) cat("\n")
    parms <- list(weights=weights, lower=lower, eps.seq=eps.seq,
                  tau=tau, Th=Th, x.center=x.center,
                  x.scale=x.scale, y.center=y.center, y.scale=y.scale,
                  monotone=monotone, additive=additive)
    parms
}
