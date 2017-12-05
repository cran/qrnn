softplus.prime <-
function(x, alpha=2)
{
    y <- exp(alpha*x)/(1+exp(alpha*x))
    ifelse(is.nan(y), 1, y)
}
