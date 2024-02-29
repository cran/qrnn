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
    ifelse(x >= 0, 1, elu(x, alpha) + alpha)
}

elu <-
function(x, alpha=1)
{
    ifelse(x >= 0, x, alpha*(exp(x)-1))
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

relu.prime <- function(x){
    ifelse(x >= 0, 1, 0)
}

relu <- function(x){ 
    ifelse(x >= 0, x, 0)
}

logistic.prime <- function(x){
    0.25/(cosh(x/2)^2)
}

logistic <- function(x){
    0.5 + 0.5*tanh(x/2)
}

lrelu.prime <- function(x){
    ifelse(x >= 0, 1, 0.01)
}

lrelu <- function(x){
    ifelse(x >= 0, x, 0.01*x)
}

softmax <- function(x){
    exp(x)/sum(exp(x))
}
