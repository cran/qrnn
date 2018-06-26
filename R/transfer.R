linear.prime <-
function(x)
{
    x*0 + 1
}

linear <-
function(x)
{
    x
}

elu.prime <-
function(x, alpha=1)
{
    y <- ifelse(x > 0, 1, elu(x, alpha) + alpha)
    ifelse(is.nan(y), 1, y)
}
elu <-
function(x, alpha=1)
{
    y <- ifelse(x > 0, x, alpha*(exp(x)-1))
    ifelse(is.infinite(y), x, y)
}
sigmoid.prime <-
function(x)
{
    (0.5)*(1-tanh(0.5*x)^2)
}

sigmoid <-
function(x)
{
    tanh(0.5*x)
}

softplus.prime <-
function(x, alpha=2)
{
    y <- exp(alpha*x)/(1+exp(alpha*x))
    ifelse(is.nan(y), 1, y)
}
softplus <-
function(x, alpha=2)
{
    y <- log(1+exp(alpha*x))/alpha
    ifelse(is.infinite(y), x, y)
}
