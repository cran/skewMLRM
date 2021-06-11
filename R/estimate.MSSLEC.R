estimate.MSSLEC <-
function(y,X=NULL,max.iter=1000,prec=1e-4,est.var=TRUE)
{
y<-as.matrix(y)
if(!is.matrix(y))
        stop("y must have at least one element")
  if(is.null(X)){X<-array(c(diag(ncol(y))),c(ncol(y),ncol(y),nrow(y)))}
  if(is.array(X)==FALSE & is.list(X)==FALSE)
        stop("X must be an array or a list")
 if(is.array(X))
  {Xs<-list()
if(ncol(y)>1 | !is.matrix(X)){
     for (i in 1:nrow(y)){
    Xs[[i]]<- matrix(t(X[,,i]),nrow=ncol(y))}}
if(ncol(y)==1 & is.matrix(X)){
     for (i in 1:nrow(y)){
    Xs[[i]]<- matrix(t(X[i,]),nrow=1)}} 
X<-Xs}
  if (ncol(y) != nrow(X[[1]]))
        stop("y does not have the same number of columns than X")
  if (nrow(y) != length(X))
        stop("y does not have the same number of observations than X")
EM.SMSNC<-function(y,X,beta0,Sigma0,eta0,nu0,dist="SSLEC",max.iter=1000,prec=prec, est.var=TRUE)
{ 
dtgamma1<-function(x,shape,rate,b=1,log=TRUE)
{
ll<-dgamma(x,shape=shape, rate=rate,log=TRUE)-pgamma(b,shape=shape,rate=rate,log.p=TRUE)
if(log==FALSE)
{ll<-exp(ll)}
ll
}
lsmsnc.prof<-function(nu,y,X,beta,Sigma,eta,dist="SSLEC")
{
aux.profile.llike.MSMSNC<-function(x,p,d,A,dist="SSLEC",nu=1){
logf<-dgamma(x,shape=(nu+p/2),rate=(d/2),log=TRUE)-pgamma(1,shape=(nu+p/2),rate=(d/2),log.p=TRUE)
log.int<-pcauchy(sqrt(x/d)*A,log.p=TRUE)+logf
exp(log.int)} 
 n=nrow(y)
  p=ncol(y)  
  m=ncol(X[[1]])
eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
ll<-c()
iSigma<-solve(Sigma)
logdSigma<-log(det(Sigma))
for(i in 1:n)
{
mu<-X[[i]]%*%beta
d=c(t(y[i,]-X[[i]]%*%beta)%*%iSigma%*%(y[i,]-X[[i]]%*%beta))
A=c(t(eta)%*%(y[i,]-X[[i]]%*%beta))
min.int<-0
max.int<-1
aux.int<-integrate(aux.profile.llike.MSMSNC,lower=min.int,upper=max.int,p=p,d=d,A=A,dist=dist,nu=nu,abs.tol=1e-15)$value
ll[i]<-log(2)+log(nu)+lgamma(nu+p/2)-(p/2)*log(2*pi)-0.5*logdSigma+pgamma(1, shape=(nu+p/2),rate=d/2,log.p=TRUE)-(nu+p/2)*log(d/2)+log(aux.int)
}
-sum(ll)
}
E.step.MSMSNC.par<-function(y,X,beta,Sigma,eta,dist="SSLEC",nu=1)
{
compute.a1<-function(a,d,p,dist="SSLEC",nu=1)
{
aux.a1.MSMSNC<-function(x,a,d,p,dist="SSLEC",nu=1){
logf<-dtgamma1(x,shape=(2*nu+p)/2,rate=d/2,log=FALSE)
log.int<-x*pcauchy(sqrt(x)*a)*logf
return(log.int)}
max.int<-1
integrate(aux.a1.MSMSNC,lower=0,upper=max.int,a,d,p,dist=dist,nu=nu,abs.tol=1e-12)$value
}
compute.b1<-function(a,d,p,dist="SSLEC",nu=1)
{
aux.b1.MSMSNC<-function(x,a,d,p,dist="SSLEC",nu=1){
logf<-dtgamma1(x,shape=(2*nu+p)/2,rate=d/2, log=FALSE)
log.int<-x*pt(sqrt(3*x)*a,df=3)*logf
return(log.int)}
max.int<-1
integrate(aux.b1.MSMSNC,lower=0,upper=max.int,a,d,p,dist=dist,nu=nu,abs.tol=1e-12)$value
}
compute.c1<-function(a,d,p,dist="SSLEC",nu=1)
{
aux.c1.MSMSNC<-function(x,a,d,p,dist="SSLEC",nu=1){
logf<-dtgamma1(x,shape=(2*nu+p)/2,rate=d/2, log=FALSE)
log.int<-sqrt(x)*dcauchy(sqrt(x)*a)*logf
return(log.int)}
max.int<-1
integrate(aux.c1.MSMSNC,lower=0,upper=max.int,a,d,p,dist=dist,nu=nu,abs.tol=1e-12)$value
}

compute.den.MSMSNC<-function(a,d,p,dist="SSLEC",nu=1)
{
 aux.den.MSMSNC<-function(x,a,d,p,dist="SSLEC",nu=1){  
logf<-dtgamma1(x,shape=(2*nu+p)/2,rate=d/2,log=FALSE)
log.int<-pcauchy(sqrt(x)*a)*logf
return(log.int)}
max.int<-1
  integrate(aux.den.MSMSNC,lower=0,upper=max.int,a,d,p,dist=dist,nu=nu,abs.tol=1e-12)$value
}
n<-nrow(y) ; p=ncol(y)
      iSigma<-solve(Sigma)
a.theta<-c();b.theta<-c();c.theta<-c()
for(i in 1:n)
{ 
      resi=y[i,]-X[[i]]%*%beta
      ai= as.numeric(t(eta)%*%resi)
      di= as.numeric(t(resi)%*%iSigma%*%resi)  
      a.a<-compute.a1(a=ai,d=di,p=p,dist=dist,nu=nu)
      a.b<-compute.b1(a=ai,d=di,p=p,dist=dist,nu=nu)
      a.c<-compute.c1(a=ai,d=di,p=p,dist=dist,nu=nu)
      den<-compute.den.MSMSNC(a=ai,d=di,p=p,dist=dist,nu=nu)
a.theta[i]<-a.a/den
b.theta[i]<-a.b/den
c.theta[i]<-as.numeric(t(eta)%*%resi)*b.theta[i]+a.c/den
}
return(list(a.theta=a.theta,b.theta=b.theta,c.theta=c.theta))
}
M1.step.MSMSNC<-function(y,X,Sigma,eta,a.theta,b.theta,c.theta)
{
p=ncol(X[[1]]);m=nrow(X[[1]])
eta<-matrix(eta,ncol=1)
sum1<-matrix(0,ncol=p,nrow=p)
sum2<-matrix(0,nrow=p,ncol=1)
  invSigma=solve(Sigma)
for(i in 1:nrow(y))
{
sum1<-sum1+t(X[[i]])%*%(a.theta[i]*invSigma+b.theta[i]*eta%*%t(eta))%*%X[[i]]
sum2<-sum2+t(X[[i]])%*%(b.theta[i]*eta%*%t(eta)%*%y[i,]+a.theta[i]*invSigma%*%y[i,] - c.theta[i]*eta)   
}
  solve(sum1,tol=1e-100)%*%sum2
}
M2.step.MSMSNC<-function(y,X,beta,a.theta)
{
p=ncol(X[[1]]);m=nrow(X[[1]]);n=nrow(y)
beta<-matrix(beta,ncol=1)
sum1<-matrix(0,ncol=m,nrow=m)
for(i in 1:n)
{
sum1<-sum1+a.theta[i]*(y[i,]-X[[i]]%*%beta)%*%t(y[i,]-X[[i]]%*%beta)
}
sum1/n
}
M3.step.MSMSNC<-function(y,X,beta,b.theta,c.theta)
{
p=ncol(X[[1]]);m=nrow(X[[1]]);n=nrow(y)
sum1<-matrix(0,ncol=m,nrow=m)
sum2<-matrix(0,ncol=1,nrow=m)
for(i in 1:n)
{
sum1<-sum1+b.theta[i]*(y[i,]-X[[i]]%*%beta)%*%t(y[i,]-X[[i]]%*%beta)
sum2<-sum2+c.theta[i]*(y[i,]-X[[i]]%*%beta)
}
solve(sum1,tol=1e-100)%*%sum2
}
n=nrow(y)
  p=ncol(y)
  m=ncol(X[[1]])
aa=system.time({
beta.last<-matrix(beta0,ncol=1)
Sigma.last<-Sigma0
eta.last<-matrix(eta0,ncol=1)
nu.last<-nu0
  lower1<-1.0001
i=0;dif=10
while(i<=max.iter & dif>prec)
{  aux<-E.step.MSMSNC.par(y,X,c(beta.last),Sigma.last,c(eta.last),dist,nu=nu.last)
   a.theta<-aux$a.theta
b.theta<-aux$b.theta
c.theta<-aux$c.theta
    beta.new<-M1.step.MSMSNC(y,X,Sigma.last,c(eta.last),a.theta,b.theta,c.theta)
   Sigma.new<-M2.step.MSMSNC(y,X,beta.new,a.theta)
eta.new<-M3.step.MSMSNC(y,X,beta.new,b.theta,c.theta)
    lambda.new<-matrix.sqrt(Sigma.new)%*%eta.new
nu.new<-optim(nu.last,lsmsnc.prof,method="Brent",lower=lower1,upper=1000,y=y,X=X,beta=beta.new,Sigma=Sigma.new,eta=eta.new,dist=dist,control = list(maxit = 1000))$par
dif=abs(lsmsnc.prof(nu.new,y,X,beta.new,Sigma.new,eta.new,dist)-lsmsnc.prof(nu.last,y,X,beta.last,Sigma.last,eta.last,dist))
eta.last=eta.new;beta.last=beta.new
Sigma.last=Sigma.new;nu.last=nu.new
i=i+1
}
 })
conv<-ifelse(i<=max.iter & dif<=prec, 0, 1)
 tempo=as.numeric(aa[3])
 lognu=-lsmsnc.prof(nu.new,y,X,beta.new,Sigma.new,eta.new,dist)
 q0=ncol(Sigma.new)
 npar=length(beta.new)+length(eta.new)+q0*(q0+1)/2+length(nu.new)
 AIC=-2*lognu+2*npar
 BIC=-2*lognu+log(nrow(y))*npar
 aux=as.list(sapply(1:p,seq,by=1,to=p))
 P=c(as.vector(beta.new),vech(matrix.sqrt(Sigma.new)),as.vector(lambda.new))
 indices=c()
 for(j in 1:p)
 {indices=c(indices,paste(j,aux[[j]],sep=""))}
 P<-matrix(P,ncol=1)
 colnames(P)<-c("estimate") 
 conv.problem=0
 se=c()
 if(est.var)
 {
 MI.obs<-FI.MSSLEC(P,y,X,nu.new)
 test=try(solve(MI.obs,tol=1e-100),silent=TRUE)
 if(is.numeric(test) & max(diag(test))<0) 
 {
 se=sqrt(-diag(test))
 }
 else conv.problem=1
 }
 nu<-nu.new
 if(est.var){P<-cbind(P,se);colnames(P)<-c("estimate","s.e.")}
rownames(P)<-c(paste("beta",1:m,sep=""),paste("alpha",indices,sep=""),paste("lambda",1:p,sep=""))
 ll<-list(estimate=P,nu=nu,logLik=lognu,AIC=AIC,BIC=BIC,iterations=i,time=tempo,conv=conv)
 if(conv.problem==1) ll$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
 ll
}
  n=nrow(y)
  p=ncol(y)
  q=ncol(X[[1]])
  b0<-matrix(0,q,q)
  b1<-matrix(0,q,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,] 
  }
  beta0<-solve(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%beta0
  }
  Sigma0<-cov(e)
  eta0<-as.matrix(moments::skewness(e))
  nu0<-3
EM.SMSNC(y,X,beta0,Sigma0,eta0,nu0,dist="SSLEC",max.iter=max.iter,prec=prec,est.var=est.var)
}
