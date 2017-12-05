qrnn.cost <-
function(weights, x, y, n.hidden, w, tau, lower, monotone, additive, eps,
         Th, Th.prime, penalty, unpenalized)
{
    penalty2 <- ifelse(identical(Th, linear), penalty, 0)
    w1w2 <- qrnn.reshape(x, y, weights, n.hidden)
    W1 <- w1w2$W1; rW1 <- nrow(W1); cW1 <- ncol(W1)
    W2 <- w1w2$W2; rW2 <- nrow(W2); cW2 <- ncol(W2)
    if (!is.null(monotone)) {
        W1[monotone,] <- exp(W1[monotone,])
        W2[1:(rW2-1), ] <- exp(W2[1:(rW2-1),])
    }
    if(!is.logical(additive)){
        W1 <- W1*additive
    }
    # Forward pass
    x <- cbind(x, 1)
    h1 <- x %*% W1
    y1 <- Th(h1)
    aug.y1 <- cbind(y1, 1)
    h2 <- aug.y1 %*% W2
    y2 <- hramp(h2, lower, eps)
    E <- y-y2
    # Backward pass
    delta2 <- hramp.prime(h2, lower, eps)*tilted.approx.prime(E, tau, eps)
    gradient.W2 <- -(t(aug.y1) %*% sweep(delta2, 1, w, '*'))
    if (!is.null(monotone)){
        gradient.W2[1:(rW2-1),] <- gradient.W2[1:(rW2-1),]*W2[1:(rW2-1),]
    }
    gradient.W2.penalty <- 2*penalty2*rbind(W2[1:(rW2-1),,drop=FALSE], 0)/
                              (length(W2)-cW2)
    E1 <- delta2 %*% t(W2[1:(rW2-1),,drop=FALSE])
    delta1 <- Th.prime(h1)*E1
    gradient.W1 = -(t(x) %*% sweep(delta1, 1, w, '*'))
    if (!is.null(monotone)){
        gradient.W1[monotone,] <- gradient.W1[monotone,]*W1[monotone,]
    }
    W1p <- W1; W1p[c(unpenalized, rW1),] <- 0
    gradient.W1.penalty <- 2*penalty*W1p/sum(W1p != 0)
    # Error & gradient
    cost <- sum(w*tilted.approx(E, tau, eps)) +
        penalty*sum(W1p^2)/sum(W1p != 0) +
        penalty2*sum(W2[1:(rW2-1),,drop=FALSE]^2)/(length(W2)-cW2)
    gradient <- c(gradient.W1 + gradient.W1.penalty,
                  gradient.W2 + gradient.W2.penalty)
    attr(cost, "gradient") <- gradient
    cost
}
