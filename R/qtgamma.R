qtgamma <-
function (p, shape, scale = 1, a = 0, b = Inf) 
{
    stopifnot(all(p >= 0 & p <= 1) & all(scale > 0) & all(shape > 
        0))
    Fa <- pgamma(a, shape, scale = scale)
    Fb <- pgamma(b, shape, scale = scale)
    pNew <- p * (Fb - Fa) + Fa
    x <- qgamma(pNew, shape, scale = scale)
    return(x)
}
