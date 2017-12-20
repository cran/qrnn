mcqrnn.fit <- function(x, y, n.hidden=2, n.hidden2=NULL, w=NULL,
                       tau=c(0.1, 0.5, 0.9), iter.max=5000,
                       n.trials=5, lower=-Inf,
                       init.range = c(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
                       monotone=NULL, eps.seq=2^(-8:-32), Th=sigmoid,
                       Th.prime=sigmoid.prime, penalty=0,
                       n.errors.max=10, trace=TRUE, ...){
    if(length(tau)==1) stop("Multiple values of \"tau\" required")
    x.y.tau <- composite.stack(x, y, tau)
    taus <- x.y.tau$tau
    xs <- cbind(taus, x.y.tau$x)
    ys <- x.y.tau$y
    if(is.null(monotone)){
        monotone <- 1
    } else{
        monotone <- c(1, monotone+1)
    }
    if(is.null(n.hidden2)){
        if(length(init.range) > 4) init.range <- init.range[1:4]
        parms <- qrnn.fit(x=xs, y=ys, n.hidden=n.hidden, w=w, tau=taus,
                          n.ensemble=1, iter.max=iter.max,
                          n.trials=n.trials, bag=FALSE, lower=lower,
                          init.range=init.range, monotone=monotone,
                          eps.seq=eps.seq, Th=Th, Th.prime=Th.prime,
                          penalty=penalty, unpenalized=1,
                          n.errors.max=n.errors.max, trace=trace, ...)
    } else{
        parms <- qrnn2.fit(x=xs, y=ys, n.hidden=n.hidden, n.hidden2=n.hidden2,
                           w=w, tau=taus, n.ensemble=1, iter.max=iter.max,
                           n.trials=n.trials, bag=FALSE, lower=lower,
                           init.range=init.range, monotone=monotone,
                           eps.seq=eps.seq, Th=Th, Th.prime=Th.prime,
                           penalty=penalty, unpenalized=1,
                           n.errors.max=n.errors.max, trace=trace, ...)
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

