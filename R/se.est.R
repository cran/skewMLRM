se.est <-
function(P,y,X,dist="MN",nu=0.5, gamma=0.5)
{
  if (!any(dist == c("MN","MT","MSL","MCN","MSN","MSNC","MSTEC","MSTN","MSTT",
"MSSLEC","MSSL","MSSL2","MSCN","MSCN2","MSCEC")))
     stop("distribution is not recognized")
y<-as.matrix(y)
if(!is.matrix(y))
        stop("y must have at least one element")
  if(is.null(X)){X<-array(c(diag(ncol(y))),c(ncol(y),ncol(y),nrow(y)))}
  if(is.array(X)==FALSE & is.list(X)==FALSE)
        stop("X must be an array or a list")
  if(is.array(X))
  {Xs<-list()
     for (i in 1:nrow(y)){
    Xs[[i]]<- matrix(t(X[,,i]),nrow=ncol(y))};X<-Xs} 
  if (ncol(y) != nrow(X[[1]]))
        stop("y does not have the same number of columns than X")
  if (nrow(y) != length(X))
        stop("y does not have the same number of observations than X")
  FI.dist <- get(paste("FI.", dist, sep = ""), mode = "function")
  P<-matrix(P,ncol=1)
if (!any(dist == c("MSTEC","MSSLEC","MSCEC"))) MI.obs<-FI.dist(P, y, X)
if (dist == "MSTEC" | dist=="MSSLEC") MI.obs<-FI.dist(P, y, X, nu)
if (dist == "MSCEC") MI.obs<-FI.dist(P, y, X, nu, gamma)
 test=try(solve(MI.obs,tol=1e-100),silent=TRUE)
 se=c()
 if(is.numeric(test) & max(diag(test))<0) 
 {
 se=sqrt(-diag(test))
 }
 else  stop("Standard errors can't be estimated: Numerical problems with the inversion of the information matrix")
 se
}
