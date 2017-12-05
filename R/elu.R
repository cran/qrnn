elu <-
function(x, alpha=1)
{
    y <- ifelse(x > 0, x, alpha*(exp(x)-1))
    ifelse(is.infinite(y), x, y)
}
