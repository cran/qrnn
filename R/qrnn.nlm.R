qrnn.nlm <-
function(x, y, n.hidden, w, tau, iter.max, n.trials, bag, lower, init.range,
         monotone, additive, eps.seq, Th, Th.prime, penalty, unpenalized,
         trace, ...)
{
    cases <- seq(nrow(x))
    if (bag) cases <- sample(nrow(x), replace=TRUE)
    x <- x[cases,,drop=FALSE]
    y <- y[cases,,drop=FALSE]
    w <- w[cases]
    if(length(tau) > 1) tau <- tau[cases]
    if(length(lower) > 1) lower <- lower[cases]
	eps.seq <- sort(eps.seq, decreasing=TRUE)
    cost.best <- Inf
    for(i in seq(n.trials)){
        weights <- qrnn.initialize(x, y, n.hidden, init.range)
        if(any(lower != -Inf)){
            for(eps in eps.seq){
                fit <- suppressWarnings(nlm(qrnn.cost, weights,
                           iterlim=iter.max, x=x, y=y, n.hidden=n.hidden,
                           w=w, tau=tau, lower=-Inf, monotone=monotone,
                           additive=additive, eps=eps, Th=Th,
                           Th.prime=Th.prime, penalty=penalty,
                           unpenalized=unpenalized,
                           check.analyticals=FALSE, ...))
                weights <- fit$estimate
            }
        } 
        for(eps in eps.seq){
            fit <- suppressWarnings(nlm(qrnn.cost, weights, iterlim=iter.max,
                       x=x, y=y, n.hidden=n.hidden, w=w, tau=tau, lower=lower,
                       monotone=monotone, additive=additive, eps=eps, Th=Th,
                       Th.prime=Th.prime, penalty=penalty,
                       unpenalized=unpenalized,
                       check.analyticals=FALSE, ...))
            weights <- fit$estimate
        }
        cost <- fit$minimum
        if(trace) cat(i, cost, "\n")
        if(cost < cost.best){
            cost.best <- cost
            weights.best <- fit$estimate
        }
    }
    if(trace) cat("*", cost.best, "\n")
    weights.best <- qrnn.reshape(x, y, weights.best, n.hidden)
    if(!is.logical(additive)){
        weights.best$W1 <- weights.best$W1*additive
    }
    weights.best
}
