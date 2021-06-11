estimate.MSSL2 <-
function(y,X=NULL,max.iter=1000,prec=1e-4,est.var=TRUE){
logsslnu<-function(nu,Y,X,beta0,Sigma,lambda){
n=nrow(Y)
p=ncol(Y)
B=matrix.sqrt(Sigma) 
Binv=solve(B)
cnf=(2*pi)^(p/2)*sqrt(det(Sigma))
auxA=t(lambda)%*%solve(B)
fy0=rep(0,n)
for(i in 1:n){
    resi<-matrix(Y[i,])-X[[i]]%*%beta0
    di<-as.numeric(t(resi)%*%Binv%*%Binv%*%resi)
    Ai=as.numeric(auxA%*%resi)
    faux <- function(u) u^(nu + p/2-1)*exp(-u*di/2)*pnorm(u^(1/2)*Ai)
    fy0[i] <- integrate(faux, 0, 1)$value}
fy=2*nu*fy0/cnf
return(sum(log(fy)))}
mssl.logL <- function(theta,Y,X){
 n=nrow(Y)
 p=ncol(Y)
 q0=ncol(X[[1]])
 pth=length(theta)-1  # especifico
 beta0=as.matrix(theta[1:q0],q0,1)
 p1=p*(p+1)/2
 B= xpnd(theta[(q0+1):(q0+p1)]) # Sigma^{1/2}
 Sigma=B%*%B
 invSigma=solve(Sigma)
 Binv=matrix.sqrt(invSigma) # Sigma^{-1/2} = B^{-1}
 lambda=as.matrix(theta[(q0+p1+1):pth],p,1)
 nu=as.numeric(theta[pth+1])
 cnf=(2*pi)^(p/2)*sqrt(det(Sigma))
 auxA=t(lambda)%*%solve(B)
 fy0=rep(0,n)
 for(i in 1:n){
    resi<-matrix(Y[i,])-X[[i]]%*%beta0
    di<-as.numeric(t(resi)%*%invSigma%*%resi)
    Ai=as.numeric(auxA%*%resi)
    faux <- function(u) u^(nu + p/2-1)*exp(-u*di/2)*pnorm(u^(1/2)*Ai)
    fy0[i] <- integrate(faux, 0, 1)$value}
fy=2*nu*fy0/cnf
return(sum(log(fy)))}
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
  beta0<-solve(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%beta0
  }
  S<-cov(e)
  invS<-solve(S)
  B<-matrix.sqrt(S)
  lambda<-as.matrix(moments::skewness(e))
  nu=5
  theta0<-c(as.vector(beta0),vech(B),as.vector(lambda),nu)
  log0= mssl.logL(theta0,y,X)
  Sigma=B%*%B
  invSigma=solve(Sigma)
  delta=lambda/sqrt(1+sum(as.vector(lambda)^2))
  Delta= matrix.sqrt(solve(Sigma))%*%delta
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
    Binv=solve(B)    
    for (i in 1:n){
          mu[i,]<-X[[i]]%*%beta0
          res[i,]<-y[i,]-X[[i]]%*%beta0
          d[i]<-as.numeric(t(res[i,])%*%Binv%*%Binv%*%res[i,])
    }
    ### Inicio passo E
    Gamainv=solve(Gama)
    A=as.vector(t(lambda)%*%solve(B)%*%t(res)) 
    
    u <- eta1 <- c()
    for (i in 1:n) {
        fauxu <- function(u) dgamma(u,shape=(2*nu+p+2)/2,rate=d[i]/2)*pnorm(u^(1/2)*A[i])
        auxu <- integrate(fauxu, 0, 1)$value/pgamma(1,shape=(2*nu+p+2)/2,rate=d[i]/2)
        faux <- function(u) u^(nu + p/2-1)*exp(-u*d[i]/2)*pnorm(u^(1/2)*A[i])
        aux22 <- integrate(faux, 0, 1)$value
        f0_fy=1/2*gamma(nu + p/2)*pgamma(1,shape=nu + p/2, rate=d[i]/2)/(((d[i]/2)^(nu + p/2))*aux22)
        u[i] = f0_fy*4*gamma((2*nu+p+2)/2)/gamma((2*nu+p)/2)/d[i]*pgamma(1,shape=(2*nu+p+2)/2,rate=d[i]/2)/pgamma(1,shape=(2*nu+p)/2,rate=d[i]/2)*auxu        
        eta1[i] =2^(nu+p/2)*gamma((2*nu+p+1)/2)*pgamma(1,shape=(2*nu+p+1)/2, rate=(d[i]+A[i]^2)/2)/(sqrt(pi)*((d[i]+A[i]^2)^((2*nu+p+1)/2))*aux22)
                }
    MT2=1/(1+as.numeric(t(Delta)%*%Gamainv%*%Delta))
    muT=MT2*as.vector(t(Delta)%*%solve(Gama)%*%t(res))
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
    beta0=solve(beta01)%*%beta02
    Delta=Delta1/sum(ut2)
    Gama=Gama0/n
    Sigma=Gama+Delta%*%t(Delta)
    invSigma=solve(Sigma)
    B=matrix.sqrt(Sigma)
    Binv=solve(B)
    delta=Binv%*%Delta
    lambda=delta/sqrt(1-as.numeric(t(delta)%*%delta))    
    #nu=optimize(logsslnu,c(1.01,40),y,X,beta0,Sigma,lambda,maximum=T) #nu,Y,X,beta0,Sigma,lambda
    #nu=as.numeric(nu$maximum)
    V<-optim(nu,logsslnu,gr=NULL,y,X,beta0,Sigma,lambda,method="L-BFGS-B",lower=1.01,upper=40,control=list(fnscale=-1,maxit=50))   
    nu<-as.numeric(V$par) 
    theta=c(as.vector(beta0),vech(B),as.vector(lambda),nu)
    logL=mssl.logL(theta,y,X)
    criterio=abs(logL-log0)
    bb=sum(eigen(Gama)$values>0)-p
    #delta=lambda/sqrt(1+as.numeric(t(lambda)%*%lambda))
    theta0=theta
    log0=logL
 }
})
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
  nu=as.numeric(theta[length(theta)])
 rownames(P)<-c(paste("beta",1:q,sep=""),paste("alpha",indices,sep=""),paste("lambda",1:p,sep=""))
 colnames(P)<-c("estimate") 
 conv.problem=0
 if(est.var)
 {
 MI.obs<-  FI.MSSL2(P,y,X,nu)
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
  object.out<-list(estimate=P,nu=nu,logLik=logL,AIC=AIC,BIC=BIC,iterations=cont,time=tempo,conv=conv)
 if(conv.problem==1) object.out$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
 object.out
 }
