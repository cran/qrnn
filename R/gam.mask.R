gam.mask <-
function(x, n.hidden)
{
    mask <- matrix(1, ncol(x)+1, n.hidden)
    nh <- ncol(mask)/ncol(x)
    ij <- matrix(seq(nh*ncol(x)), nrow=nh)
    for(i in seq(ncol(x))){
        mask[i,c(ij[,-i])] <- 0
    }
    mask
}
