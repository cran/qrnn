huber.prime <-
function(x, eps)
{
    dh <- x/eps
    dh[x>eps] <- 1
    dh[x< -eps] <- -1
    dh[is.nan(dh)] <- 0
    dh
}

huber <-
function(x, eps)
{
    h <- ifelse(abs(x)>eps, abs(x)-eps/2, (x^2)/(2*eps))
    h[is.nan(h)] <- 0
    h
}

tilted.approx.prime <-
function(x, tau, eps)
{
    ifelse(x>0, tau*huber.prime(x, eps), (1-tau)*huber.prime(x, eps))
}
tilted.approx <-
function(x, tau, eps)
{
    ifelse(x>0, tau*huber(x, eps), (1-tau)*huber(x, eps))
}
hramp.prime <-
function(x, lower, eps)
{
    if(length(lower) > 1){
        mapply(hramp.prime, x, lower, eps)
    } else{
        if (lower==-Inf){
            return(1)
        } else{
            dhr <- (x-lower)/eps
            dhr[x>(lower+eps)] <- 1
            dhr[x<lower] <- 0
            return(dhr)
        }
    }
}
hramp <-
function(x, lower, eps)
{
    if(length(lower) > 1){
        mapply(hramp, x, lower, eps)
    } else{
        if (lower==-Inf){
            return(x)
        } else{
            return(ifelse(x>lower, huber(x-lower, eps), 0)+lower)
        }
    }
}
