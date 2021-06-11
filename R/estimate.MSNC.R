estimate.MSNC <-
function(y,X=NULL,max.iter=1000,prec=1e-4,est.var=TRUE)
{
lmsncr<- function(y,X=NULL,beta0,Sigma,lambda)
{
  if(is.array(X)==FALSE & is.list(X)==FALSE)
        stop("X must be an array or a list")
  y<-as.matrix(y)
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
        stop("y does not have the observations than X")  
n=nrow(y)
  eta=solve(matrix.sqrt(Sigma))%*%lambda
  logv=0
  for (i in 1:n){
      mui<-X[[i]]%*%beta0
      ai= as.numeric(t(eta)%*%(y[i,]-mui))
      aux<-2*mvtnorm::dmvnorm(y[i,],as.vector(mui),sigma=Sigma)*pcauchy(ai)
      logv=logv+log(aux)
    } 
 return(logv)
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
aa=system.time({
  m=ncol(X[[1]])
  p=nrow(X[[1]])
  n=length(X)
  b0<-matrix(0,m,m)
  b1<-matrix(0,m,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,]
  }
  beta.last<-solve(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%beta.last
  }
  Sigma.last<-cov(e)
  Sm.inv<-solve(matrix.sqrt(Sigma.last))
  lambda.last<-as.matrix(moments::skewness(e))
  eta.last=Sm.inv%*%lambda.last
  
  log0<- lmsncr(y,X,beta.last,Sigma.last,lambda.last)

iter=0;dif=1
while(iter<=max.iter && dif>prec)
{
a.theta<-rep(1,n)                               
b.theta<-c()                           
c.theta<-c() 
    for (i in 1:n){
         ai= as.numeric(t(eta.last)%*%(y[i,]-X[[i]]%*%beta.last))
         b.theta[i]<-pt(sqrt(3)*ai,df=3)/pcauchy(ai)
         c.theta[i]<-ai*b.theta[i]+dcauchy(ai)/pcauchy(ai)
    }                        
beta.new<-M1.step.MSMSNC(y,X,Sigma.last,c(eta.last),a.theta,b.theta,c.theta) 
Sigma.new<-M2.step.MSMSNC(y,X,beta.new,a.theta)  
eta.new<-M3.step.MSMSNC(y,X,beta.new,b.theta,c.theta) 
    lambda.new<-matrix.sqrt(Sigma.new)%*%eta.new
    lognu<-lmsncr(y,X,beta.new,Sigma.new,lambda.new)
    dif<-abs(lognu-log0)
eta.last=eta.new;beta.last=beta.new;Sigma.last=Sigma.new
iter=iter+1
    log0=lognu 
}
 })
 tempo=as.numeric(aa[3])
 npar=length(beta.new)+length(eta.new)+p*(p+1)/2
 AIC=-2*lognu+2*npar
 BIC=-2*lognu+log(n)*npar
 P<-matrix(c(as.vector(beta.new),vech(matrix.sqrt(Sigma.new)),as.vector(lambda.new)),ncol=1)
 colnames(P)<-c("estimate") 
 conv.problem=0
 if(est.var)
 {
 MI.obs<-FI.MSNC(P,y,X)
 test=try(solve(MI.obs,tol=1e-100),silent=TRUE)
 se=c()
 if(is.numeric(test) & max(diag(test))<0) 
 {
 se=sqrt(-diag(test))
 P<-cbind(P,se)
 colnames(P)<-c("estimate","s.e.")
 }
 else conv.problem=1
 }
conv<-ifelse(iter<=max.iter & dif<=prec, 0, 1)
 aux=as.list(sapply(1:p,seq,by=1,to=p))
 indices=c()
 for(j in 1:p)
 {indices=c(indices,paste(j,aux[[j]],sep=""))}
 rownames(P)<-c(paste("beta",1:m,sep=""),paste("alpha",indices,sep=""),paste("lambda",1:p,sep=""))
 ll<-list(estimate=P,logLik=lognu,AIC=AIC,BIC=BIC,iterations=iter,time=tempo, conv=conv)
 if(conv.problem==1) ll$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
 object.out<-ll
object.out
}
