rMSL <-
function(n, mu,Sigma, nu=1){
   if(is.null(n))
    stop("sample size must be specified")
  if(is.null(mu))
    stop("mu must be specified")
  if(is.null(Sigma))
    stop("Sigma must be specified")
    if(is.vector(Sigma)) Sigma=matrix(Sigma)
  if(round(n)!=n | n<=0)
    stop("sample size must be a positive integer")      

  if(!is.matrix(mu)) mu<-matrix(mu, ncol=length(mu),nrow=1)
  if (!isSymmetric(Sigma))
        stop("Sigma is not symmetric")
  if (min(eigen(Sigma)$values)<=0)
        stop("Sigma is not a positive-definite matrix")
  if (ncol(Sigma) != ncol(mu))
        stop("The dimension of mu does not agree with the dimension of Sigma")
  if(nu<=0)
        stop("nu must be positive")
  p<-ncol(Sigma)
  y<-matrix(0,n,p)
 for(i in 1:n){
           u<-rbeta(1,nu,1)
      y[i,]=as.vector(mu)+MASS::mvrnorm(1,rep(0,p),Sigma/u)
  } 
    return(y)
  }
