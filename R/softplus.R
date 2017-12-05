softplus <-
function(x, alpha=2)
{
    y <- log(1+exp(alpha*x))/alpha
    ifelse(is.infinite(y), x, y)
}
