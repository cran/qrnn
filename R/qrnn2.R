qrnn2.cost <-
function(weights, x, y, n.hidden, n.hidden2, w, tau, lower, monotone, eps,
         Th, Th.prime, penalty, unpenalized)
{
    w1w2w3 <- qrnn2.reshape(x, y, weights, n.hidden, n.hidden2)
    W1 <- w1w2w3$W1; rW1 <- nrow(W1); cW1 <- ncol(W1)
    W2 <- w1w2w3$W2; rW2 <- nrow(W2); cW2 <- ncol(W2)
    W3 <- w1w2w3$W3; rW3 <- nrow(W3); cW3 <- ncol(W3)
    if (!is.null(monotone)) {
        W1[monotone,] <- exp(W1[monotone,])
        W2[1:(rW2-1), ] <- exp(W2[1:(rW2-1),])
        W3[1:(rW3-1), ] <- exp(W3[1:(rW3-1),])
    }
    # Forward pass
    x <- cbind(x, 1)
    h1 <- x %*% W1
    y1 <- Th(h1)
    aug.y1 <- cbind(y1, 1)
    h2 <- aug.y1 %*% W2
    y2 <- Th(h2)
    aug.y2 <- cbind(y2, 1)
    h3 <- aug.y2 %*% W3
    y3 <- hramp(h3, lower, eps)
    E <- y-y3
    # Backward pass
    delta3 <- hramp.prime(h3, lower, eps)*tilted.approx.prime(E, tau, eps)
    gradient.W3 <- -(t(aug.y2) %*% sweep(delta3, 1, w, "*"))
    if (!is.null(monotone)){
        gradient.W3[1:(rW3-1),] <- gradient.W3[1:(rW3-1),]*W3[1:(rW3-1),]
    }
    E2 <- delta3 %*% t(W3[1:(rW3-1),,drop=FALSE])
    delta2 <- Th.prime(h2)*E2
    gradient.W2 <- -(t(aug.y1) %*% sweep(delta2, 1, w, "*"))
    if (!is.null(monotone)){
        gradient.W2[1:(rW2-1),] <- gradient.W2[1:(rW2-1),]*W2[1:(rW2-1),]
    }
    gradient.W2.penalty <- 2*penalty*rbind(W2[1:(rW2-1),,drop=FALSE], 0)/
                              (length(W2)-cW2)
    E1 <- delta2 %*% t(W2[1:(rW2-1),,drop=FALSE])
    delta1 <- Th.prime(h1)*E1
    gradient.W1 = -(t(x) %*% sweep(delta1, 1, w, "*"))
    if (!is.null(monotone)){
        gradient.W1[monotone,] <- gradient.W1[monotone,]*W1[monotone,]
    }
    W1p <- W1; W1p[c(unpenalized, rW1),] <- 0
    gradient.W1.penalty <- 2*penalty*W1p/sum(W1p != 0)
    # Error & gradient
    cost <- sum(w*tilted.approx(E, tau, eps)) +
        penalty*sum(W1p^2)/sum(W1p != 0) +
        penalty*sum(W2[1:(rW2-1),,drop=FALSE]^2)/(length(W2)-cW2)
    gradient <- c(gradient.W1 + gradient.W1.penalty,
                  gradient.W2 + gradient.W2.penalty,
                  gradient.W3)
    attr(cost, "gradient") <- gradient
    cost
}
qrnn2.eval <-
function(x, W1, W2, W3, lower, monotone, eps, Th)
{
    if (!is.null(monotone)) {
        W1[monotone, ] <- exp(W1[monotone, ])
        W2[1:(nrow(W2) - 1), ] <- exp(W2[1:(nrow(W2) - 1),])
        W3[1:(nrow(W3) - 1), ] <- exp(W3[1:(nrow(W3) - 1),])
    }
    x <- cbind(x, 1)
    h1 <- x %*% W1
    y1 <- Th(h1)
    aug.y1 <- cbind(y1, 1)
    h2 <- aug.y1 %*% W2
    y2 <- Th(h2)
    aug.y2 <- cbind(y2, 1)
    h3 <- aug.y2 %*% W3
    y3 <- hramp(h3, lower, eps)
    y3
}
qrnn2.initialize <-
function(x, y, n.hidden, n.hidden2,
         init.range=c(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5))
{
    if(length(init.range)==6){
        r11 <- init.range[1]
        r12 <- init.range[2]
        r21 <- init.range[3]
        r22 <- init.range[4]
        r31 <- init.range[5]
        r32 <- init.range[6]
    } else{
        r11 <- r21 <- r31 <- init.range[1]
        r12 <- r22 <- r32 <- init.range[2]
    }
    W1 <- matrix(runif((ncol(x)+1)*n.hidden, r11, r12), ncol(x)+1, n.hidden)
    W2 <- matrix(runif((n.hidden+1)*n.hidden2, r21, r22), n.hidden+1, n.hidden2)
    W3 <- matrix(runif((n.hidden2+1)*ncol(y), r31, r32), n.hidden2+1, ncol(y))
    c(W1, W2, W3)
}
adam <- function(f, p, x, y, w, tau, ..., iterlim=5000, iterbreak=iterlim,
                 alpha=0.01, minibatch=nrow(x), beta1=0.9, beta2=0.999,
                 epsilon=1e-8, print.level=10){
    minibatches <- suppressWarnings(matrix(seq_along(y), nrow=minibatch))
    f.best <- f(p, x=x, y=y, w=w, tau=tau, ...)
    p.best <- p
    i.break <- 0
    M <- R <- p*0
    if(print.level > 0) cat(0, f.best, i.break, f.best, "\n")
    for(iter in seq(iterlim)){
        cases.random <- sample(nrow(x))
        for(i in seq(ncol(minibatches))){
            cases <- cases.random[minibatches[,i]]
            grad <- attr(f(p, x=x[cases,,drop=FALSE], y=y[cases,,drop=FALSE],
                         w=w[cases], tau=tau[cases], ...), "gradient")
            M <- beta1*M + (1-beta1)*grad
            R <- beta2*R + (1-beta2)*grad^2
            m_k_hat <- M/(1-beta1^iter)
            r_k_hat <- R/(1-beta2^iter)
            p <- p - alpha*m_k_hat/(sqrt(r_k_hat) + epsilon)
        }
        f.iter <- f(p, x=x, y=y, w=w, tau=tau, ...)
        if(f.iter < f.best){
            f.best <- f.iter
            p.best <- p
            i.break <- 0
        } else{
            i.break <- i.break + 1
        }
        if(print.level > 0 && iter%%print.level==0)
            cat(iter, f.iter, i.break, f.best, "\n")
        if(i.break > iterbreak){
            break
        }
    }
    if(print.level > 0) cat("****", f.best, "\n")
    list(estimate=p.best, minimum=f.best)
}
qrnn2.fit <-
function(x, y, n.hidden=2, n.hidden2=2, w=NULL, tau=0.5, n.ensemble=1,
         iter.max=5000, n.trials=5, bag=FALSE, lower=-Inf,
         init.range=c(-0.5, 0.5, -0.5, 0.5, -0.5, 0.5),
         monotone=NULL, eps.seq=2^(-8:-32), Th=sigmoid, Th.prime=sigmoid.prime,
         penalty=0, unpenalized=NULL, n.errors.max=10, trace=TRUE,
         method=c("nlm", "adam"), ...)

{
    method <- match.arg(method)
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    if (!is.matrix(y)) stop("\"y\" must be a matrix")
    if (ncol(y) != 1) stop("\"y\" must be univariate")
    if (any(is.na(c(x, y)))) stop("missing values in \"x\" or \"y\"")
    if (any(apply(x, 2, sd) < .Machine$double.eps^0.5)) stop("zero variance column(s) in \"x\"")
    if (any((tau > 1) | (tau < 0))) stop("invalid \"tau\"")
    if (identical(Th, linear))
        stop("use \"qrnn.fit\" for linear models")
    is.whole <- function(x, tol = .Machine$double.eps^0.5) 
        abs(x - round(x)) < tol
    if(is.null(w)) w <- rep(1/nrow(y), nrow(y))
    if (any(w < 0)) stop("invalid \"w\"")
    x <- scale(x)
    x.center <- attr(x, "scaled:center")
    x.scale <- attr(x, "scaled:scale")
    y <- scale(y)
    y.center <- attr(y, "scaled:center")
    y.scale <- attr(y, "scaled:scale")
    lower.scaled <- (lower-y.center)/y.scale
    weights <- vector("list", n.ensemble)
    if(trace) cat("tau =", unique(tau), "\n", sep=" ")
    for (i in seq(n.ensemble)){
        if(trace) cat(i, "/", n.ensemble, "\n", sep="")
        w.tmp <- NA
        class(w.tmp) <- "try-error"
        n.errors <- 0
        while (inherits(w.tmp, "try-error")) {
            w.tmp <- try(qrnn2.optimize(x, y, n.hidden, n.hidden2, w, tau,
                             iter.max, n.trials, bag, lower.scaled,
                             init.range, monotone, eps.seq, Th, Th.prime,
                             penalty, unpenalized, trace, method, ...),
                        silent = TRUE)
            n.errors <- n.errors + 1
            if (n.errors > n.errors.max) stop("optimization failed")
        }
        weights[[i]] <- w.tmp
    }
    if(trace) cat("\n")
    parms <- list(weights=weights, lower=lower, eps.seq=eps.seq,
                  tau=tau, Th=Th, x.center=x.center,
                  x.scale=x.scale, y.center=y.center, y.scale=y.scale,
                  monotone=monotone)
    parms
}
qrnn2.optimize <-
function(x, y, n.hidden, n.hidden2, w, tau, iter.max, n.trials, bag, lower,
         init.range, monotone, eps.seq, Th, Th.prime, penalty, unpenalized,
         trace, method, ...)
{
    cases <- seq(nrow(x))
    if (bag) cases <- sample(nrow(x), replace=TRUE)
    x <- x[cases,,drop=FALSE]
    y <- y[cases,,drop=FALSE]
    w <- w[cases]
    if(length(lower) > 1) lower <- lower[cases]
    if(length(tau)==1) tau <- rep(tau, length(y))
    tau <- tau[cases]
	eps.seq <- sort(eps.seq, decreasing=TRUE)
    cost.best <- Inf
    for(i in seq(n.trials)){
        weights <- qrnn2.initialize(x, y, n.hidden, n.hidden2, init.range)
        if(any(lower != -Inf)){
            for(eps in eps.seq){
                if(method=="nlm"){
                    fit <- suppressWarnings(nlm(qrnn2.cost, weights,
                               iterlim=iter.max, x=x, y=y, n.hidden=n.hidden,
                               n.hidden2=n.hidden2, w=w, tau=tau, lower=-Inf,
                               monotone=monotone, eps=eps, Th=Th,
                               Th.prime=Th.prime, penalty=penalty,
                               unpenalized=unpenalized,
                               check.analyticals=FALSE, ...))
                } else if(method=="adam"){
                    fit <- suppressWarnings(adam(qrnn2.cost, weights, x, y,
                               w, tau, iterlim=iter.max, n.hidden=n.hidden,
                               n.hidden2=n.hidden2, lower=-Inf,
                               monotone=monotone, eps=eps, Th=Th,
                               Th.prime=Th.prime, penalty=penalty,
                               unpenalized=unpenalized, ...))
                }
                weights <- fit$estimate
            }
        } 
        for(eps in eps.seq){
            if(method=="nlm"){
                fit <- suppressWarnings(nlm(qrnn2.cost, weights,
                           iterlim=iter.max, x=x, y=y, n.hidden=n.hidden,
                           n.hidden2=n.hidden2, w=w, tau=tau, lower=lower,
                           monotone=monotone, eps=eps, Th=Th, Th.prime=Th.prime,
                           penalty=penalty, unpenalized=unpenalized,
                           check.analyticals=FALSE, ...))
            } else if(method=="adam"){
                fit <- suppressWarnings(adam(qrnn2.cost, weights, x, y, w, tau,
                           iterlim=iter.max, n.hidden=n.hidden,
                           n.hidden2=n.hidden2, lower=lower, monotone=monotone,
                           eps=eps, Th=Th, Th.prime=Th.prime,
                           penalty=penalty, unpenalized=unpenalized,
                           minibatch=nrow(x), ...))
            }
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
    weights.best <- qrnn2.reshape(x, y, weights.best, n.hidden, n.hidden2)
    weights.best
}
qrnn2.predict <-
function(x, parms)
{
    if (!is.matrix(x)) stop("\"x\" must be a matrix")
    weights <- parms$weights
    lower <- parms$lower
    monotone <- parms$monotone
    eps <- min(parms$eps.seq)
    Th <- parms$Th
    x.center <- parms$x.center
    x.scale <- parms$x.scale
    y.center <- parms$y.center
    y.scale <- parms$y.scale
    lower <- (lower-y.center)/y.scale
    x <- sweep(x, 2, x.center, "-")
    x <- sweep(x, 2, x.scale, "/")
    y.bag <- matrix(0, ncol=length(weights), nrow=nrow(x))
    for (i in seq_along(weights)){
        y.bag[,i] <- qrnn2.eval(x, weights[[i]]$W1, weights[[i]]$W2,
                               weights[[i]]$W3, lower, monotone,
                               eps, Th)
        y.bag[,i] <- y.bag[,i]*y.scale + y.center
    }
    y.bag
}
qrnn2.reshape <-
function(x, y, weights, n.hidden, n.hidden2)
{
    N11 <- ncol(x)+1
    N12 <- n.hidden
    N1 <- N11*N12
    W1 <- weights[1:N1]
    W1 <- matrix(W1, N11, N12)
    N21 <- n.hidden+1
    N22 <- n.hidden2
    N2 <- N1 + N21*N22
    W2 <- weights[(N1+1):N2]
    W2 <- matrix(W2, N21, N22)
    N31 <- n.hidden2+1
    N32 <- ncol(y)
    N3 <- N2 + N31*N32
    W3 <- weights[(N2+1):N3]
    W3 <- matrix(W3, N31, N32)
    list(W1=W1, W2=W2, W3=W3)
}
