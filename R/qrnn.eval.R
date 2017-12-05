qrnn.eval <-
function(x, W1, W2, lower, monotone, additive, eps, Th)
{
    if (!is.null(monotone)) {
        W1[monotone, ] <- exp(W1[monotone, ])
        W2[1:(nrow(W2) - 1), ] <- exp(W2[1:(nrow(W2) - 1),])
    }
    if(!is.logical(additive)){
        W1 <- W1*additive
    }
    x <- cbind(x, 1)
    h1 <- x %*% W1
    y1 <- Th(h1)
    aug.y1 <- cbind(y1, 1)
    y2 <- aug.y1 %*% W2
    y2 <- hramp(y2, lower, eps)
    y2
}
