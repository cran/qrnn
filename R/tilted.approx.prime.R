tilted.approx.prime <-
function(x, tau, eps)
{
    ifelse(x>0, tau*huber.prime(x, eps), (1-tau)*huber.prime(x, eps))
}
