FI.MSCEC <-
function(P,y,X,nu,gamma)
{
    y <- as.matrix(y)
    if (!is.matrix(y)) 
        stop("y must have at least one element")
        if (missing(X)) {
            X <- array(1, c(ncol(y), 1, nrow(y)))
        }
    if (is.null(X)) {
        X <- array(c(diag(ncol(y))), c(ncol(y), ncol(y), nrow(y)))
    }
    if (is.array(X) == FALSE & is.list(X) == FALSE) 
        stop("X must be an array or a list")
    if (is.array(X)) {
        Xs <- list()
        for (i in 1:nrow(y)) {
            Xs[[i]] <- matrix(t(X[, , i]), nrow = ncol(y))
        }
        X <- Xs
    }
    if (ncol(y) != nrow(X[[1]])) 
        stop("y does not have the same number of columns than X")
    if (nrow(y) != length(X)) 
        stop("y does not have the same number of observations than X")
    if (nu < 0 | nu > 1) 
        stop("nu must be between 0 and 1")
    if (gamma < 0 | gamma > 1) 
        stop("gamma must be between 0 and 1")
FI.MSMSNC<- function(P,y,X,nu,dist=c("STEC","SSLEC","SCEC"))
{
dtgamma1<-function(x,shape,rate,b=1,log=TRUE)
{
ll<-dgamma(x,shape=shape, rate=rate,log=TRUE)-pgamma(b,shape=shape,rate=rate,log.p=TRUE)
if(log==FALSE)
{ll<-exp(ll)}
ll
}
Iphi.MSMSNC<-function(a,d,p,w,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
aux.Iphi.MSMSNC<-function(x,a,d,p,w,dist=c("STEC","SSLEC","SCEC"),nu=1){
if(dist=="STEC" || dist=="SSLEC")
logf<-switch(dist,
STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2),
SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2,log=FALSE))
log.int<-x^w*exp(-0.5*x*d)*dcauchy(sqrt(x)*a)*logf
return(log.int)}
    if(dist=="STEC" | dist=="SSLEC"){
max.int<-switch(dist,STEC=Inf,SSLEC=1)}
integrate(aux.Iphi.MSMSNC,lower=0,upper=max.int,a,d,p,w,dist=dist,nu=nu,abs.tol=1e-12)$value
}
IPhi.MSMSNC<-function(a,d,p,w,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
aux.IPhi.MSMSNC<-function(x,a,d,p,w,dist=c("STEC","SSLEC","SCEC"),nu=1){
if(dist=="STEC" || dist=="SSLEC")
logf<-switch(dist,
STEC=dgamma(x,shape=(nu+p)/2,rate=(nu+d)/2),
SSLEC=dtgamma1(x,shape=(2*nu+p)/2,rate=d/2,log=FALSE))
log.int<-x^w*exp(-0.5*x*d)*pcauchy(sqrt(x)*a)*logf
return(log.int)}
if(dist=="STEC" | dist=="SSLEC") max.int<-switch(dist,STEC=Inf,SSLEC=1)
integrate(aux.IPhi.MSMSNC,lower=0,upper=max.int,a,d,p,w,dist=dist,nu=nu,abs.tol=1e-12)$value
}
 if(is.array(X))
  {Xs<-list()
if(ncol(y)>1 | !is.matrix(X)){
     for (i in 1:nrow(y)){
    Xs[[i]]<- matrix(t(X[,,i]),nrow=ncol(y))}}
if(ncol(y)==1 & is.matrix(X)){
     for (i in 1:nrow(y)){
    Xs[[i]]<- matrix(t(X[i,]),nrow=1)}} 
X<-Xs}
 # theta=list(beta=beta.new,Sigma=Sigma.new,eta=eta.new,nu=nu,dif=dif,iter=i)
  p=ncol(y)
  n=nrow(y)
  m=ncol(X[[1]])
  beta=P[1:m]
  q0=length(beta)
  Sigma<-matrix(0,ncol=p, nrow=p)
  Sigma[lower.tri(Sigma, diag = TRUE)]<-P[1:(p*(p+1)/2)+m]
  Sigma[upper.tri(Sigma, diag = TRUE)]<-P[1:(p*(p+1)/2)+m]
  Sinv=solve2(Sigma) # Sigma^{-1}
  C=matrix.sqrt(Sigma)
  Cinv=solve2(C)# Sigma^{-1/2}
  eta=P[m+(p*(p+1)/2)+1:p]
  lambda=C%*%eta
  #Delta=Sminv%*%lambda  # eta
  q1=p*(p+1)/2
  pth=q0+p+q1
  dlogS=matrix(0,pth,1)
        # dlog Sigma
      for (r in 1:q1){
        indr=rep(0,q1)
        indr[r]=1
        Cr=xpnd(indr)
        auxi=Sinv%*%(Cr%*%C+C%*%Cr)
        dlogS[q0+r]=sum(diag(auxi))
        }

  dAi=matrix(0,pth,1)
  ddi=matrix(0,pth,1)
  MI=matrix(0,pth,pth)
  for (i in 1:n){
      resi=y[i,]-X[[i]]%*%beta  # p x 1
      ai= as.numeric(t(eta)%*%resi)
      di= as.numeric(t(resi)%*%Sinv%*%resi)
      if(dist!="SCEC")
    {
      Iphi.MSMSNC.i<-Iphi.MSMSNC(a=ai,d=di,p=p,w=(p+1)/2,dist=dist,nu=nu)
      IPhi.MSMSNC.i<-IPhi.MSMSNC(a=ai,d=di,p=p,w=(p+2)/2,dist=dist,nu=nu)
      Ki<-IPhi.MSMSNC(a=ai,d=di,p=p,w=p/2,dist=dist,nu=nu)
      }
      if(dist=="SCEC"){  
         nu1=nu[1]
         nu2=nu[2]
         Iphi.MSMSNC.i<-nu1*(nu2^((p+1)/2))*exp(-0.5*nu2^0.5*di)*dcauchy(nu2^0.25*ai)+(1-nu1)*exp(-0.5*di)*dcauchy(ai)
         IPhi.MSMSNC.i<-nu1*(nu2^((p+2)/2))*exp(-0.5*nu2^0.5*di)*pcauchy(nu2^0.25*ai)+(1-nu1)*exp(-0.5*di)*pcauchy(ai)
         Ki<-nu1*(nu2^(p/2))*exp(-0.5*nu2^0.5*di)*pcauchy(nu2^0.25*ai)+(1-nu1)*exp(-0.5*di)*pcauchy(ai)
      }
      dAi[1:q0]=t(X[[i]])%*%Cinv%*%lambda
      dAi[(q0+q1+1):pth]=Cinv%*%resi
      ddi[1:q0]=-2*t(X[[i]])%*%Cinv%*%Cinv%*%resi
      for (r in 1:q1){
           indr=rep(0,q1)
           indr[r]=1
           Cr=xpnd(indr)
           ddi[q0+r]=-t(resi)%*%Cinv%*%(Cr%*%C+C%*%Cr)%*%Cinv%*%resi
           dAi[q0+r]=-t(lambda)%*%Cinv%*%Cr%*%Cinv%*%resi
           }
      scoretheta=-0.5*dlogS +(Iphi.MSMSNC.i*dAi-0.5*IPhi.MSMSNC.i*ddi)/Ki
      MI=MI+scoretheta%*%t(scoretheta)
    } 
 return(-MI)
}
nu=c(nu,gamma)
 FI.MSMSNC(P,y,X,dist="SCEC",nu)
}
