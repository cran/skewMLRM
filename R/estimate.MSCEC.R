estimate.MSCEC <-
function(y,X=NULL,max.iter=1000,prec=1e-4,est.var=TRUE, nu.fixed=0.1, gamma.fixed=0.5)
{
y.or<-y; X.or<-X
y<-as.matrix(y)
if(nu.fixed!=FALSE & !is.numeric(nu.fixed))
   stop("nu fixed must be a number between 0 and 1")
if(is.numeric(nu.fixed) & (as.numeric(nu.fixed)<0 | as.numeric(nu.fixed)>1))
   stop("nu fixed must be a number between 0 and 1")
if(gamma.fixed!=FALSE & !is.numeric(gamma.fixed))
   stop("gamma fixed must be a number between 0 and 1")
if(is.numeric(gamma.fixed) & (as.numeric(gamma.fixed)<0 | as.numeric(gamma.fixed)>1))
   stop("gamma fixed must be a number between 0 and 1")
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
logscecgama<-function(gama,p,d,V,A){
fy=p/2*sum(V)*log(gama)-1/2*gama*sum(d*V)+sum(V*log(pcauchy(sqrt(gama)*A)))
return(fy)}
EM.SMSNC<-function(y,X,beta0,Sigma0,eta0,nu0,dist=c("STEC","SSLEC","SCEC"),max.iter=1000,prec=prec, est.var=TRUE)
{ 
dtgamma1<-function(x,shape,rate,b=1,log=TRUE)
{
ll<-dgamma(x,shape=shape, rate=rate,log=TRUE)-pgamma(b,shape=shape,rate=rate,log.p=TRUE)
if(log==FALSE)
{ll<-exp(ll)}
ll
}
lsmsnc.prof<-function(nu,y,X,beta,Sigma,eta,dist="SCEC")
{
 n=nrow(y)
  p=ncol(y)  
  m=ncol(X[[1]])
eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
nu1<-nu[1];nu2<-nu[2]
ll<-c()
Sigma2<-Sigma/nu2
for(i in 1:n)
{
mu<-X[[i]]%*%beta
A=c(t(eta)%*%(y[i,]-X[[i]]%*%beta))
ll[i]<-log(2*nu1*mvtnorm::dmvnorm(y[i,],mu,Sigma2)*pcauchy(sqrt(nu2)*A)+2*(1-nu1)*mvtnorm::dmvnorm(y[i,],mu,Sigma)*pcauchy(A))
}
-sum(ll)
}
E.step.MSMSNC.par<-function(y,X,beta,Sigma,eta,dist=c("STEC","SSLEC","SCEC"),nu=1)
{
compute.abc.SCEC<-function(y,X,beta,Sigma,eta,nu)
{ 
  m=length(y) 
  p=length(beta)
nu1<-nu[1];nu2<-nu[2]
eta<-matrix(eta,ncol=1);beta<-matrix(beta,ncol=1)
y<-matrix(y,ncol=1)
Sigma2<-Sigma/nu2
mu<-c(X%*%beta)
  aux <-(nu1*mvtnorm::dmvnorm(c(y),mu,Sigma2)+(1-nu1)*mvtnorm::dmvnorm(c(y),mu,Sigma))
  q<-nu1*mvtnorm::dmvnorm(c(y),mu,Sigma2)/(nu1*mvtnorm::dmvnorm(c(y),mu,Sigma2)+(1-nu1)*mvtnorm::dmvnorm(c(y),mu,Sigma))
  res=y-X%*%beta
  den<-pcauchy(sqrt(nu2)*c(t(eta)%*%res))*q+pcauchy(c(t(eta)%*%res))*(1-q)
num.a<-nu2*pcauchy(sqrt(nu2)*c(t(eta)%*%res))*q+pcauchy(c(t(eta)%*%res))*(1-q)
num.b<-nu2*pt(sqrt(3*nu2)*c(t(eta)%*%res),df=3)*q+pt(sqrt(3)*c(t(eta)%*%res),df=3)*(1-q)
num.c<-sqrt(nu2)*dcauchy(sqrt(nu2)*c(t(eta)%*%res))*q+dcauchy(c(t(eta)%*%res))*(1-q)
a<-num.a/den
b<-num.b/den
c0<-t(eta)%*%res*b+num.c/den
return(list(a.theta=a,b.theta=b,c.theta=c(c0)))
}
n<-nrow(y) ; p=ncol(y)
      a.theta<-c();b.theta<-c();c.theta<-c()
for(i in 1:n)
{
aux<-compute.abc.SCEC(y[i,],X[[i]],beta,Sigma,eta,nu)
a.theta[i]<-aux$a.theta
b.theta[i]<-aux$b.theta
c.theta[i]<-aux$c.theta
}
return(list(a.theta=a.theta,b.theta=b.theta,c.theta=c.theta))
}
M1.step.MSMSNC<-function(y,X,Sigma,eta,a.theta,b.theta,c.theta)
{
p=ncol(X[[1]]);m=nrow(X[[1]])
eta<-matrix(eta,ncol=1)
sum1<-matrix(0,ncol=p,nrow=p)
sum2<-matrix(0,nrow=p,ncol=1)
  invSigma=solve2(Sigma)
for(i in 1:nrow(y))
{
sum1<-sum1+t(X[[i]])%*%(a.theta[i]*invSigma+b.theta[i]*eta%*%t(eta))%*%X[[i]]
sum2<-sum2+t(X[[i]])%*%(b.theta[i]*eta%*%t(eta)%*%y[i,]+a.theta[i]*invSigma%*%y[i,] - c.theta[i]*eta)   
}
  solve2(sum1)%*%sum2
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
solve2(sum1)%*%sum2
}
n=nrow(y)
  p=ncol(y)
  m=ncol(X[[1]])
if(!nu.fixed & !gamma.fixed){
aa=system.time({
beta.last<-matrix(beta0,ncol=1)
Sigma.last<-Sigma0
eta.last<-matrix(eta0,ncol=1)
  lambda.last <- matrix.sqrt(Sigma.last)%*%eta.last
nu.last<-nu0
  nu=nu.last[1]
  gama=nu.last[2]
  mu=matrix(0,n,p)
  d<-rep(0,n)
  res<-matrix(0,n,p)
  lower1<-switch(dist, STEC=2.001, SSLEC=1.001)
i=0;dif=10
while(i<=max.iter && dif>prec)
{ aux<-E.step.MSMSNC.par(y,X,c(beta.last),Sigma.last,c(eta.last),dist,nu=nu.last)
   a.theta<-aux$a.theta
b.theta<-aux$b.theta
c.theta<-aux$c.theta
   invSigma=solve2(Sigma.last)
   for (i in 1:n){
          mu[i,]<-X[[i]]%*%beta.last
          res[i,]<-y[i,]-X[[i]]%*%beta.last
          d[i]<-as.numeric(t(res[i,])%*%invSigma%*%res[i,])
    }
    B=matrix.sqrt(Sigma.last)
    A=as.vector(t(lambda.last)%*%solve2(B)%*%t(res))
    fy=2*nu*mvtnorm::dmvnorm(res, sigma = Sigma.last/gama)*pcauchy(sqrt(gama)*A)+2*(1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma.last)*pcauchy(A)
    V= 2*nu*mvtnorm::dmvnorm(res, sigma = Sigma.last/gama)*pcauchy(sqrt(gama)*A)/fy
    beta.new<-M1.step.MSMSNC(y,X,Sigma.last,c(eta.last),a.theta,b.theta,c.theta)
   Sigma.new<-M2.step.MSMSNC(y,X,beta.new,a.theta)
eta.new<-M3.step.MSMSNC(y,X,beta.new,b.theta,c.theta)
    lambda.new<-matrix.sqrt(Sigma.new)%*%eta.new
    nu=mean(V)
gamaf<-optim(gama,logscecgama,gr = NULL,method="Brent",lower=0.001,upper=0.999,p,d,V,A,control = list(maxit = 1000,fnscale=-1))$par
      gama=as.numeric(gamaf) 
      nu.new=c(nu,gama)
dif=abs(lsmsnc.prof(nu.new,y,X,beta.new,Sigma.new,eta.new,dist)-lsmsnc.prof(nu.last,y,X,beta.last,Sigma.last,eta.last,dist))
eta.last=eta.new;beta.last=beta.new;lambda.last=lambda.new
Sigma.last=Sigma.new;nu.last=nu.new
i=i+1
}
 })}
if(!nu.fixed & is.numeric(gamma.fixed)){
aa=system.time({
beta.last<-matrix(beta0,ncol=1)
Sigma.last<-Sigma0
eta.last<-matrix(eta0,ncol=1)
  lambda.last <- matrix.sqrt(Sigma.last)%*%eta.last
nu.last<-nu0
  nu=nu.last[1]
  gama=gamma.fixed
  mu=matrix(0,n,p)
  d<-rep(0,n)
  res<-matrix(0,n,p)
  lower1<-switch(dist, STEC=2.001, SSLEC=1.001)
i=0;dif=10
while(i<=max.iter && dif>prec)
{ aux<-E.step.MSMSNC.par(y,X,c(beta.last),Sigma.last,c(eta.last),dist,nu=nu.last)
   a.theta<-aux$a.theta
b.theta<-aux$b.theta
c.theta<-aux$c.theta
   invSigma=solve2(Sigma.last)
   for (i in 1:n){
          mu[i,]<-X[[i]]%*%beta.last
          res[i,]<-y[i,]-X[[i]]%*%beta.last
          d[i]<-as.numeric(t(res[i,])%*%invSigma%*%res[i,])
    }
    B=matrix.sqrt(Sigma.last)
    A=as.vector(t(lambda.last)%*%solve2(B)%*%t(res))
    fy=2*nu*mvtnorm::dmvnorm(res, sigma = Sigma.last/gama)*pcauchy(sqrt(gama)*A)+2*(1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma.last)*pcauchy(A)
    V= 2*nu*mvtnorm::dmvnorm(res, sigma = Sigma.last/gama)*pcauchy(sqrt(gama)*A)/fy
    beta.new<-M1.step.MSMSNC(y,X,Sigma.last,c(eta.last),a.theta,b.theta,c.theta)
   Sigma.new<-M2.step.MSMSNC(y,X,beta.new,a.theta)
eta.new<-M3.step.MSMSNC(y,X,beta.new,b.theta,c.theta)
    lambda.new<-matrix.sqrt(Sigma.new)%*%eta.new
    nu=mean(V)
      nu.new=c(nu,gama)
dif=abs(lsmsnc.prof(nu.new,y,X,beta.new,Sigma.new,eta.new,dist)-lsmsnc.prof(nu.last,y,X,beta.last,Sigma.last,eta.last,dist))
eta.last=eta.new;beta.last=beta.new;lambda.last=lambda.new
Sigma.last=Sigma.new;nu.last=nu.new
i=i+1
}
 })}
if(is.numeric(nu.fixed) & !gamma.fixed){
aa=system.time({
beta.last<-matrix(beta0,ncol=1)
Sigma.last<-Sigma0
eta.last<-matrix(eta0,ncol=1)
  lambda.last <- matrix.sqrt(Sigma.last)%*%eta.last
nu.last<-nu0
  nu=nu.fixed
  gama=nu.last[2]
  mu=matrix(0,n,p)
  d<-rep(0,n)
  res<-matrix(0,n,p)
  lower1<-switch(dist, STEC=2.001, SSLEC=1.001)
i=0;dif=10
while(i<=max.iter && dif>prec)
{ aux<-E.step.MSMSNC.par(y,X,c(beta.last),Sigma.last,c(eta.last),dist,nu=nu.last)
   a.theta<-aux$a.theta
b.theta<-aux$b.theta
c.theta<-aux$c.theta
   invSigma=solve2(Sigma.last)
   for (i in 1:n){
          mu[i,]<-X[[i]]%*%beta.last
          res[i,]<-y[i,]-X[[i]]%*%beta.last
          d[i]<-as.numeric(t(res[i,])%*%invSigma%*%res[i,])
    }
    B=matrix.sqrt(Sigma.last)
    A=as.vector(t(lambda.last)%*%solve2(B)%*%t(res))
    fy=2*nu*mvtnorm::dmvnorm(res, sigma = Sigma.last/gama)*pcauchy(sqrt(gama)*A)+2*(1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma.last)*pcauchy(A)
    V= 2*nu*mvtnorm::dmvnorm(res, sigma = Sigma.last/gama)*pcauchy(sqrt(gama)*A)/fy
    beta.new<-M1.step.MSMSNC(y,X,Sigma.last,c(eta.last),a.theta,b.theta,c.theta)
   Sigma.new<-M2.step.MSMSNC(y,X,beta.new,a.theta)
eta.new<-M3.step.MSMSNC(y,X,beta.new,b.theta,c.theta)
    lambda.new<-matrix.sqrt(Sigma.new)%*%eta.new
gamaf<-optim(gama,logscecgama,gr = NULL,method="Brent",lower=0.001,upper=0.999,p,d,V,A,control = list(maxit = 1000,fnscale=-1))$par
      gama=as.numeric(gamaf) 
      nu.new=c(nu,gama)
dif=abs(lsmsnc.prof(nu.new,y,X,beta.new,Sigma.new,eta.new,dist)-lsmsnc.prof(nu.last,y,X,beta.last,Sigma.last,eta.last,dist))
eta.last=eta.new;beta.last=beta.new;lambda.last=lambda.new
Sigma.last=Sigma.new;nu.last=nu.new
i=i+1
}
 })}
if(is.numeric(nu.fixed) & is.numeric(gamma.fixed)){
aa=system.time({
beta.last<-matrix(beta0,ncol=1)
Sigma.last<-Sigma0
eta.last<-matrix(eta0,ncol=1)
  lambda.last <- matrix.sqrt(Sigma.last)%*%eta.last
nu.last<-nu0
  nu=nu.fixed
  gama=gamma.fixed
  mu=matrix(0,n,p)
  d<-rep(0,n)
  res<-matrix(0,n,p)
  lower1<-switch(dist, STEC=2.001, SSLEC=1.001)
i=0;dif=10
while(i<=max.iter && dif>prec)
{ aux<-E.step.MSMSNC.par(y,X,c(beta.last),Sigma.last,c(eta.last),dist,nu=nu.last)
   a.theta<-aux$a.theta
b.theta<-aux$b.theta
c.theta<-aux$c.theta
   invSigma=solve2(Sigma.last)
   for (i in 1:n){
          mu[i,]<-X[[i]]%*%beta.last
          res[i,]<-y[i,]-X[[i]]%*%beta.last
          d[i]<-as.numeric(t(res[i,])%*%invSigma%*%res[i,])
    }
    B=matrix.sqrt(Sigma.last)
    A=as.vector(t(lambda.last)%*%solve2(B)%*%t(res))
    fy=2*nu*mvtnorm::dmvnorm(res, sigma = Sigma.last/gama)*pcauchy(sqrt(gama)*A)+2*(1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma.last)*pcauchy(A)
    beta.new<-M1.step.MSMSNC(y,X,Sigma.last,c(eta.last),a.theta,b.theta,c.theta)
   Sigma.new<-M2.step.MSMSNC(y,X,beta.new,a.theta)
eta.new<-M3.step.MSMSNC(y,X,beta.new,b.theta,c.theta)
    lambda.new<-matrix.sqrt(Sigma.new)%*%eta.new
      nu.new=c(nu,gama)
dif=abs(lsmsnc.prof(nu.new,y,X,beta.new,Sigma.new,eta.new,dist)-lsmsnc.prof(nu.last,y,X,beta.last,Sigma.last,eta.last,dist))
eta.last=eta.new;beta.last=beta.new;lambda.last=lambda.new
Sigma.last=Sigma.new;nu.last=nu.new
i=i+1
}
 })}
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
 conv.problem=1
 if(est.var)
 {
 MI.obs<-FI.MSCEC(P,y,X,nu.new[1],nu.new[2])
 test=try(solve2(MI.obs),silent=TRUE)
 if(is.numeric(test) & max(diag(test))<0) 
 {
 conv.problem=0
 se=sqrt(-diag(test))
 }
 }
 nu<-nu.new[1]
 gamma<-nu.new[2]
 if(est.var){P<-cbind(P,se); colnames(P)<-c("estimate","s.e.")}
rownames(P)<-c(paste("beta",1:m,sep=""),paste("alpha",indices,sep=""),paste("lambda",1:p,sep=""))
 if(conv.problem==0)  ll<-list(coefficients=P[,1],se=P[,2],nu=nu,gamma=gamma,logLik=lognu,AIC=AIC,BIC=BIC,iterations=i,time=tempo,conv=conv,dist="MSCEC",class="MSMSNC",n=nrow(y))
 else{
 ll<-list(coefficients=P[,1],nu=nu,gamma=gamma,logLik=lognu,AIC=AIC,BIC=BIC,iterations=i,time=tempo,conv=conv,dist="MSCEC",class="MSMSNC",n=nrow(y))
 ll$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
 }
 class(ll) <- "skewMLRM"
 ll$y<-y.or
 ll$X<-X.or
 ll$"function"<-"estimate.MSCEC"
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
  beta0<-solve2(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%beta0
  }
  Sigma0<-cov(e)
  eta0<-as.matrix(moments::skewness(e))
  nu0<-c(0.5,0.5)
EM.SMSNC(y,X,beta0,Sigma0,eta0,nu0,dist="SCEC",max.iter=max.iter,prec=prec,est.var=est.var)
}
