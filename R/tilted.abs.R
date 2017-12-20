tilted.abs <-
function(x, tau)
{
    ifelse(x>0, x*tau, x*(tau-1))
}

