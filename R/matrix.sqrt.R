matrix.sqrt <-
function(A)
{
   sva <- svd(A)
    if (min(sva$d)>=0)
       Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
    else
       stop("Matrix square root is not defined")
    return(Asqrt)
}
