estimate.MSTT <-
function(y,X,max.iter=1000,prec=1e-4,est.var=TRUE,nu.fixed=3,nu.min=2.0001){
y.or<-y; X.or<-X
logsttnu<-function(nu,Y,X,beta0,Sigma,lambda){
n=nrow(Y)
p=ncol(Y)
B=matrix.sqrt(Sigma) 
Binv=solve2(B)
mu=matrix(0,n,p)
d<-matrix(0,n,1)
res<-matrix(0,n,p)
for(i in 1:n){
    mu[i,]<-X[[i]]%*%beta0
    res[i,]<-Y[i,]-X[[i]]%*%beta0
    d[i]<-t(res[i,])%*%Binv%*%Binv%*%res[i,]
}
A=as.vector(t(lambda)%*%solve2(B)%*%t(res))
fy=2*mvtnorm::dmvt(res, sigma = Sigma, df = nu, log = FALSE)*pt(sqrt((nu+p)/(nu+d))*A, nu+p)
return(sum(log(fy)))}
mstt.logL <- function(theta,Y,X){
 n=nrow(Y)
 p=ncol(Y)
 q0=ncol(X[[1]])
 pth=length(theta)-1  # especifico
 beta0=as.matrix(theta[1:q0],q0,1)
 p1=p*(p+1)/2
 B= xpnd(theta[(q0+1):(q0+p1)]) # Sigma^{1/2}
 Sigma=B%*%B
 invSigma=solve2(Sigma)
 Binv=matrix.sqrt(invSigma) # Sigma^{-1/2} = B^{-1}
 lambda=as.matrix(theta[(q0+p1+1):pth],p,1)
 nu=as.numeric(theta[pth+1])
 mu=matrix(0,n,p)
 d<-matrix(0,n,1)
 res<-matrix(0,n,p)
 for(i in 1:n){
    mu[i,]<-X[[i]]%*%beta0
    res[i,]<-Y[i,]-X[[i]]%*%beta0
    d[i]<-t(res[i,])%*%Binv%*%Binv%*%res[i,]}
A=as.vector(t(lambda)%*%solve2(B)%*%t(res))
fy=2*mvtnorm::dmvt(res, sigma = Sigma, df = nu, log = FALSE)*pt(sqrt((nu+p)/(nu+d))*A, nu+p) 
return(sum(log(fy)))
}
y<-as.matrix(y)
if(!is.numeric(nu.min) | nu.min<=0) stop("nu.min should be a positive number")
if(nu.fixed!=FALSE & !is.numeric(nu.fixed))
   stop("nu fixed must be a number greater than 1")
if(is.numeric(nu.fixed) & as.numeric(nu.fixed)<1)
   stop("nu fixed must be a number greater than 1")
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
if(!nu.fixed){
  aa=system.time({
  n=nrow(y)
  p=ncol(y)
  #array(1,c(p,1,n))
  q=ncol(X[[n]])
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
  S<-cov(e)
  invS<-solve2(S)
  B<-matrix.sqrt(S)
  lambda<-as.matrix(moments::skewness(e))
  nu=5
  theta0<-c(as.vector(beta0),vech(B),as.vector(lambda),nu)
  log0= mstt.logL(theta0,y,X)
  Sigma=B%*%B
  #SinvSigma=solve2(Sigma)
  delta=lambda/sqrt(1+sum(as.vector(lambda)^2))
  Delta= matrix.sqrt(solve2(Sigma))%*%delta
  a1=Sigma-Delta%*%t(Delta)
  bb=eigen(a1)$values
  { if (sum(bb>0)==p) Gama=a1
  else Gama=Sigma
  }
  criterio=1
  cont=0
  bb=1
  mu=matrix(0,n,p)
  d<-rep(0,n)
  res<-matrix(0,n,p)
  while((criterio>=prec)&&(cont<=max.iter)&&(bb>=0)){
    cont=cont+1
    Binv=solve2(B)    
    for (i in 1:n){
          mu[i,]<-X[[i]]%*%beta0
          res[i,]<-y[i,]-X[[i]]%*%beta0
          d[i]<-as.numeric(t(res[i,])%*%Binv%*%Binv%*%res[i,])
    }
    ### Inicio passo E
    Gamainv=solve2(Gama)
    A=as.vector(t(lambda)%*%solve2(B)%*%t(res))    
    u=2*gamma((nu+p+2)/2)/gamma((nu+p)/2)/(nu+d)*pt(sqrt((nu+p+2)/(nu+d))*A, nu+p+2)/pt(sqrt((nu+p)/(nu+d))*A, nu+p)
    eta1=1/pt(sqrt((nu+p)/(nu+d))*A, nu+p)/sqrt(pi)*gamma((nu+p+1)/2)/gamma((nu+p)/2)*((nu+d)^((nu+p)/2))/((nu+d+A^2)^((nu+p+1)/2))
    MT2=1/(1+as.numeric(t(Delta)%*%Gamainv%*%Delta))
    muT=MT2*as.vector(t(Delta)%*%solve2(Gama)%*%t(res))
    ut = u*muT+sqrt(MT2)*eta1
    ut2=u*(muT^2)+MT2+muT*sqrt(MT2)*eta1
    ### Fim passo E
    beta01=matrix(0,q,q)
    beta02=matrix(0,q,1)
    Delta1=matrix(0,p,1)    
    Gama0=matrix(0,p,p)
    for (i in 1:n){
        yi=as.matrix(y[i,])
        Xi=X[[i]]  # p x q
        resi=as.matrix(res[i,])        
        beta01=beta01+u[i]*t(Xi)%*%Gamainv%*%Xi
        beta02=beta02+t(Xi)%*%Gamainv%*%(u[i]*yi-ut[i]*Delta)
        Delta1=Delta1+ut[i]*resi
        Gama0=Gama0+u[i]*resi%*%t(resi)-ut[i]*(Delta%*%t(resi)+resi%*%t(Delta))+ut2[i]*Delta%*%t(Delta)
    }
    beta0=solve2(beta01)%*%beta02
    Delta=Delta1/sum(ut2)
    Gama=Gama0/n
    Sigma=Gama+Delta%*%t(Delta)
    #invSigma=solve2(Sigma)
    B=matrix.sqrt(Sigma)
    Binv=solve2(B)
    delta=Binv%*%Delta
    lambda=delta/sqrt(1-as.numeric(t(delta)%*%delta))    
    nu=optimize(logsttnu,c(nu.min,100),y,X,beta0,Sigma,lambda,maximum=T) 
    nu=as.numeric(nu$maximum) 
    theta=c(as.vector(beta0),vech(B),as.vector(lambda),nu)
    logL=mstt.logL(theta,y,X)
    criterio=abs(logL-log0)
    bb=sum(eigen(Gama)$values>0)-p
    #delta=lambda/sqrt(1+as.numeric(t(lambda)%*%lambda))
    theta0=theta
    log0=logL
 }
})}
if(is.numeric(nu.fixed)){
  aa=system.time({
  n=nrow(y)
  p=ncol(y)
  #array(1,c(p,1,n))
  q=ncol(X[[n]])
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
  S<-cov(e)
  invS<-solve2(S)
  B<-matrix.sqrt(S)
  lambda<-as.matrix(moments::skewness(e))
  nu=nu.fixed
  theta0<-c(as.vector(beta0),vech(B),as.vector(lambda),nu)
  log0= mstt.logL(theta0,y,X)
  Sigma=B%*%B
  #SinvSigma=solve2(Sigma)
  delta=lambda/sqrt(1+sum(as.vector(lambda)^2))
  Delta= matrix.sqrt(solve2(Sigma))%*%delta
  a1=Sigma-Delta%*%t(Delta)
  bb=eigen(a1)$values
  { if (sum(bb>0)==p) Gama=a1
  else Gama=Sigma
  }
  criterio=1
  cont=0
  bb=1
  mu=matrix(0,n,p)
  d<-rep(0,n)
  res<-matrix(0,n,p)
  while((criterio>=prec)&&(cont<=max.iter)&&(bb>=0)){
    cont=cont+1
    Binv=solve2(B)    
    for (i in 1:n){
          mu[i,]<-X[[i]]%*%beta0
          res[i,]<-y[i,]-X[[i]]%*%beta0
          d[i]<-as.numeric(t(res[i,])%*%Binv%*%Binv%*%res[i,])
    }
    ### Inicio passo E
    Gamainv=solve2(Gama)
    A=as.vector(t(lambda)%*%solve2(B)%*%t(res))    
    u=2*gamma((nu+p+2)/2)/gamma((nu+p)/2)/(nu+d)*pt(sqrt((nu+p+2)/(nu+d))*A, nu+p+2)/pt(sqrt((nu+p)/(nu+d))*A, nu+p)
    eta1=1/pt(sqrt((nu+p)/(nu+d))*A, nu+p)/sqrt(pi)*gamma((nu+p+1)/2)/gamma((nu+p)/2)*((nu+d)^((nu+p)/2))/((nu+d+A^2)^((nu+p+1)/2))
    MT2=1/(1+as.numeric(t(Delta)%*%Gamainv%*%Delta))
    muT=MT2*as.vector(t(Delta)%*%solve2(Gama)%*%t(res))
    ut = u*muT+sqrt(MT2)*eta1
    ut2=u*(muT^2)+MT2+muT*sqrt(MT2)*eta1
    ### Fim passo E
    beta01=matrix(0,q,q)
    beta02=matrix(0,q,1)
    Delta1=matrix(0,p,1)    
    Gama0=matrix(0,p,p)
    for (i in 1:n){
        yi=as.matrix(y[i,])
        Xi=X[[i]]  # p x q
        resi=as.matrix(res[i,])        
        beta01=beta01+u[i]*t(Xi)%*%Gamainv%*%Xi
        beta02=beta02+t(Xi)%*%Gamainv%*%(u[i]*yi-ut[i]*Delta)
        Delta1=Delta1+ut[i]*resi
        Gama0=Gama0+u[i]*resi%*%t(resi)-ut[i]*(Delta%*%t(resi)+resi%*%t(Delta))+ut2[i]*Delta%*%t(Delta)
    }
    beta0=solve2(beta01)%*%beta02
    Delta=Delta1/sum(ut2)
    Gama=Gama0/n
    Sigma=Gama+Delta%*%t(Delta)
    #invSigma=solve2(Sigma)
    B=matrix.sqrt(Sigma)
    Binv=solve2(B)
    delta=Binv%*%Delta
    lambda=delta/sqrt(1-as.numeric(t(delta)%*%delta))    
    theta=c(as.vector(beta0),vech(B),as.vector(lambda),nu)
    logL=mstt.logL(theta,y,X)
    criterio=abs(logL-log0)
    bb=sum(eigen(Gama)$values>0)-p
    #delta=lambda/sqrt(1+as.numeric(t(lambda)%*%lambda))
    theta0=theta
    log0=logL
 }
})}
  npar=length(theta)
  AIC=-2*logL+2*npar
  BIC=-2*logL+log(n)*npar
  tempo=as.numeric(aa[3])
 conv<-ifelse(cont<=max.iter & criterio<=prec, 0, 1)
   aux=as.list(sapply(1:p,seq,by=1,to=p))
 indices=c()
 for(j in 1:p)
 {indices=c(indices,paste(j,aux[[j]],sep=""))}
 P<-matrix(theta[-length(theta)],ncol=1)
 nu<-theta[length(theta)]
 rownames(P)<-c(paste("beta",1:q,sep=""),paste("alpha",indices,sep=""),paste("lambda",1:p,sep=""))
 colnames(P)<-c("estimate") 
 conv.problem=1
 if(est.var)
 {
 MI.obs<-  FI.MSTT(P,y,X,nu)
 test=try(solve2(MI.obs),silent=TRUE)
 se=c()
 if(is.numeric(test) & max(diag(test))<0) 
 {
 conv.problem=0
 se=sqrt(-diag(test))
 P<-cbind(P,se)
 colnames(P)<-c("estimate","s.e.")
 }
 }
if(conv.problem==0) ll<-list(coefficients=P[,1],se=P[,2],nu=nu,logLik=logL,AIC=AIC,BIC=BIC,iterations=cont,time=tempo,conv=conv,dist="MSTT",class="MSMSN",n=nrow(y))
else{
ll<-list(coefficients=P[,1],nu=nu,logLik=logL,AIC=AIC,BIC=BIC,iterations=cont,time=tempo,conv=conv,dist="MSTT",class="MSMSN",n=nrow(y))
 ll$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
}
 object.out<-ll
 class(object.out) <- "skewMLRM"
 object.out$y<-y.or
 object.out$X<-X.or
 object.out$"function"<-"estimate.MSTT"
object.out
 }
