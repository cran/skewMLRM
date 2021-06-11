dtgamma <-
function (x, shape, scale = 1, a = 0, b = Inf) 
{
    stopifnot(all(shape > 0) & all(scale > 0))
    Fa <- pgamma(a, shape, scale = scale)
    Fb <- pgamma(b, shape, scale = scale)
    y <- dgamma(x, shape, scale = scale)
    inda <- which(x < a)
    indb <- which(x > b)
    if (length(inda) > 0) 
        y[inda] <- 0
    if (length(indb) > 0) 
        y[indb] <- 0
    return(y/(Fb - Fa))
}
