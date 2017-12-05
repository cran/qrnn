qrnn.initialize <-
function(x, y, n.hidden, init.range=c(-0.5, 0.5, -0.5, 0.5))
{
    if(length(init.range)==4){
        r11 <- init.range[1]
        r12 <- init.range[2]
        r21 <- init.range[3]
        r22 <- init.range[4]
    } else{
        r11 <- r21 <- init.range[1]
        r12 <- r22 <- init.range[2]
    }
    W1 <- matrix(runif((ncol(x)+1)*n.hidden, r11, r12), ncol(x)+1, n.hidden)
    W2 <- matrix(runif((n.hidden+1)*ncol(y), r21, r22), n.hidden+1, ncol(y))
    c(W1, W2)
}
