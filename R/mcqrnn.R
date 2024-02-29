mcqrnn.fit <- function(x, y, n.hidden=2, n.hidden2=NULL, w=NULL,
                       tau=c(0.1, 0.5, 0.9), iter.max=5000,
                       n.trials=5, lower=-Inf,
                       init.range = c(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
                       monotone=NULL, eps.seq=2^seq(-8, -32, by=-4),
                       Th=sigmoid, Th.prime=sigmoid.prime, penalty=0,
                       n.errors.max=10, trace=TRUE,
                       method=c("nlm", "adam"), scale.y=TRUE, ...){
    if(length(tau)==1 && is.integer(tau)){
        if(trace) cat(paste("Stochastic estimation of quantile regression process using", tau, "samples\n"))
        xs <- matrix(NA, ncol=ncol(x), nrow=tau)
        ys <- matrix(NA, ncol=1, nrow=tau)
        taus <- rep(NA, tau)
        if(!is.null(w)) ws <- rep(NA, tau)
        for(i in seq(tau)){
            case.i <- sample(nrow(x), size=1)
            x.i <- x[case.i,,drop=FALSE]
            y.i <- y[case.i,,drop=FALSE]
            tau.i <- runif(1)
            x.y.tau.i <- composite.stack(x.i, y.i, tau.i)
            xs[i,] <- x.y.tau.i$x
            ys[i,] <- x.y.tau.i$y
            taus[i] <- x.y.tau.i$tau
            if(!is.null(w)) ws[i] <- w[case.i]
        }
        xs <- cbind(taus, xs)
        if(!is.null(w)) w <- ws
    } else{
        if(length(tau)==1)
            stop("Improper \'tau\' for stochastic estimation of quantile regression process (e.g., \'tau\' is not an integer)")
        x.y.tau <- composite.stack(x, y, tau)
        taus <- x.y.tau$tau
        xs <- cbind(taus, x.y.tau$x)
        ys <- x.y.tau$y
    }
    if(is.null(monotone)){
        monotone <- 1
    } else{
        monotone <- c(1, monotone+1)
    }
    if(is.null(n.hidden2)){
        if(!is.list(init.range) && length(init.range) > 4)
            init.range <- init.range[1:4]
        parms <- qrnn.fit(x=xs, y=ys, n.hidden=n.hidden, w=w, tau=taus,
                          n.ensemble=1, iter.max=iter.max,
                          n.trials=n.trials, bag=FALSE, lower=lower,
                          init.range=init.range, monotone=monotone,
                          eps.seq=eps.seq, Th=Th, Th.prime=Th.prime,
                          penalty=penalty, unpenalized=1,
                          n.errors.max=n.errors.max, trace=trace,
                          scale.y=scale.y, ...)
    } else{
        parms <- qrnn2.fit(x=xs, y=ys, n.hidden=n.hidden, n.hidden2=n.hidden2,
                           w=w, tau=taus, n.ensemble=1, iter.max=iter.max,
                           n.trials=n.trials, bag=FALSE, lower=lower,
                           init.range=init.range, monotone=monotone,
                           eps.seq=eps.seq, Th=Th, Th.prime=Th.prime,
                           penalty=penalty, unpenalized=1,
                           n.errors.max=n.errors.max, trace=trace, 
                           method=method, scale.y=scale.y, ...)
    }
    parms
}

mcqrnn.predict <- function(x, parms, tau=NULL){
    if(is.null(tau)) tau <- unique(parms$tau)
    x.tau.stack <- composite.stack(x, NA, tau)
    taus <- x.tau.stack$tau
    xs <- cbind(taus, x.tau.stack$x)
    if(any(grepl('W3', sapply(parms$weights, names)))){
        pred <- matrix(qrnn2.predict(xs, parms), ncol=length(tau))
    } else{
        pred <- matrix(qrnn.predict(xs, parms), ncol=length(tau))
    }
    colnames(pred) <- paste0('tau=', tau)
    pred
}