elu.prime <-
function(x, alpha=1)
{
    y <- ifelse(x > 0, 1, elu(x, alpha) + alpha)
    ifelse(is.nan(y), 1, y)
}
