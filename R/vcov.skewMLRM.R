vcov.skewMLRM <-
function(object, ...) 
{
if(object$"function"!="distMahal")
{
FI.dist <- get(paste("FI.", object$dist, sep = ""), mode = "function")
P<-coef(object); dist=object$dist
if (!any(dist == c("MSTEC","MSSLEC","MSCEC","MSTT","MSSL2","MSCN2"))) MI.obs<-FI.dist(P=P, y=object$y, X=object$X)
if (dist == "MSTEC" | dist=="MSSLEC" | dist=="MSTT" | dist=="MSSL2") MI.obs<-FI.dist(P=P, y=object$y, X=object$X, nu=object$nu)
if (dist == "MSCEC" | dist=="MSCN2") MI.obs<-FI.dist(P=P, y=object$y, X=object$X, nu=object$nu, gamma=object$gamma)
 test=try(solve2(MI.obs),silent=TRUE)
if(is.numeric(test) & max(diag(test))<0) 
 {
  colnames(test)=rownames(test)=names(coef(object))
 }
 else stop("Variance-covariance matrix can't be estimated: Numerical problems with the inversion of the information matrix")
}
if(object$"function"=="distMahal")
{
cat("\nvcov can't be used with an skewMLRM object created with distMahal function", sep = "")
}
-test
}
