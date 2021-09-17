solve2 <-
function(A)
{
A<-as.matrix(A)
if(!is.matrix(A)) stop("A should be a matrix")
if(ncol(A)!=nrow(A)) stop("A should be a square matrix")
A_lu <- matrixcalc::lu.decomposition(A)
backsolve(A_lu$U, backsolve(A_lu$L, diag(ncol(A)), upper.tri=F))
}
