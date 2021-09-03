estimate.MSL <-
function(y, X, max.iter=1000, prec=1e-4, est.var=TRUE, nu.min=2.0001)
{
y.or<-y;X.or<-X
fauxs<-function(u,v,p,d,r,s){u^(v+p/2-r)*exp(-u*d/2)*log(u)^s}
lmsr<-function(P,y,X){
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p*(p+1)/2)])
  invB=solve(B)
  v=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]}
  l=sum(log(v/(det(B)*(2*pi)^(p/2))*gamma(v+p/2)*pgamma(1,v+p/2,scale=2/d)/((d/2)^(v+p/2))))
  return(l)}
lmsvr<-function(v,b,B,y,X){
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  invB=solve(B)
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]}
  l=sum(log(v/(det(B)*(2*pi)^(p/2))*gamma(v+p/2)*pgamma(1,v+p/2,scale=2/d)/((d/2)^(v+p/2))))
  return(l)}  
elogu<-function(u,nu,p,d){log(u)*dtgamma(u, shape=nu+p/2, scale=2/d, a=0, b=1)}
escmsr<-function(P,y,X){
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  p0=p*(p+1)/2
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=matrix(0, nrow=p)
  v=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  W<-matrix(0,n,1)
  u<-matrix(0,n,1)
  dldb1<-matrix(0,q,1)
  dldb2<-matrix(0,q,1)
  dlda<-matrix(0,p0,1)
  dldv0<-0
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(2*v+p)/d[i]*pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i])
    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
    W[i]<-dnorm(auxi)/pnorm(auxi)
    dldb1<-dldb1+as.numeric(u[i])*t(X[[i]])%*%invB%*%invB%*%e[i,]
    dldb2<-dldb2+as.numeric(W[i])*t(X[[i]])%*%invB%*%h
    ddida<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    dldda<-matrix(0,p0,1)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      dldda[j]<-2*sum(diag(invB%*%Bj))
      dAida[j]<--t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      ddida[j]<--t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]}
    dlda<-dlda-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(W[i])*dAida
    dldv0<-dldv0+integrate(fauxs,lower=0,upper=1,v,p,d[i],1,1)$val/integrate(fauxs,lower=0,upper=1,v,p,d[i],1,0)$val}
  dldb=dldb1-dldb2
  dldv=dldv0+n/v
  dldh=invB%*%t(e)%*%W
  el=c(as.vector(dldb),as.vector(dlda),as.vector(dldh),dldv)
  return(el[-c(q+p0+1:p)])} 
y<-as.matrix(y)
if(!is.numeric(nu.min) | nu.min<=0) stop("nu.min should be a positive number")
if(!is.matrix(y)) stop("y must have at least one element")
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
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b0<-matrix(0,q,q)
  b1<-matrix(0,q,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,]
  }
  b<-solve(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
  }
  S<-cov(e)
  invS<-solve(S)
  B<-matrix.sqrt(S)
  invB<-solve(B)
  h<-matrix(0, nrow=p)
  v<-1.01
  P_0<-c(as.vector(b),vech(B),as.vector(h),v)
  log0=lmsr(P_0,y,X)
  crit<-1
  iter<-0
  while ((crit>prec)&&(iter<=max.iter))
  { 
    iter<-iter+1 
    d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b
      d[i]<-t(e[i,])%*%invS%*%e[i,]
    }
    u<-((2*v+p)/d)*(pgamma(1,v+1+p/2,scale=2/d)/pgamma(1,v+p/2,scale=2/d))
    aux<-pmax(as.vector(e%*%invB%*%h),-37)
    W<-as.matrix(dnorm(aux)/pnorm(aux))
    t<-e%*%invB%*%h+W 
    S<-1/n*(t(e)%*%diag(as.vector(u))%*%e)
    invS<-solve(S)
    B<-matrix.sqrt(S)
    invB<-solve(B)
    b0<-matrix(0,q,q)
    b1<-matrix(0,q,1)
    for(i in 1:n){
      b0<-b0+t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
      b1<-b1+t(X[[i]])%*%(as.numeric(u[i])*invS%*%y[i,]-as.numeric(t[i])*invB%*%h+invB%*%h%*%t(h)%*%invB%*%y[i,])  
    }
    b<-solve(b0)%*%b1
    V<-optim(v,lmsvr,gr=NULL,b,B,y,X,method="L-BFGS-B",lower=nu.min,upper=100,control=list(fnscale=-1,maxit=50))   
    v<-as.numeric(V$par)
    P<-c(as.vector(b),vech(B),as.vector(h),v)
    logvero=lmsr(P,y,X)
    crit=abs(logvero-log0)
    log0=logvero
  }
  })
  npar=length(P)
  AIC=-2*logvero+2*npar
  BIC=-2*logvero+log(nrow(y))*npar
 conv<-ifelse(iter<=max.iter & crit<=prec, 0, 1)
  tempo=as.numeric(aa[3])
  aux=as.list(sapply(1:p,seq,by=1,to=p))
 P<-matrix(P[-c(q+p*(p+1)/2+1:p)],ncol=1)
 colnames(P)<-c("estimate") 
 conv.problem=1
 if(est.var)
 {
 MI.obs<-FI.MSL(P,y,X)
 test=try(solve(MI.obs,tol=1e-100),silent=TRUE)
 se=c()
 if(is.numeric(test) & max(diag(test))<0) 
 {
 conv.problem=0
 se=sqrt(-diag(test))
 P<-cbind(P,se)
 colnames(P)<-c("estimate","s.e.")
 }
 }
 indices=c()
 for(j in 1:p)
 {indices=c(indices,paste(j,aux[[j]],sep=""))}
 rownames(P)<-c(paste("beta",1:q,sep=""),paste("alpha",indices,sep=""),"nu")
if(conv.problem==0)  ll<-list(coefficients=P[,1],se=P[,2],logLik=logvero,AIC=AIC,BIC=BIC,iterations=iter,time=tempo,conv=conv,dist="MSL",class="MSMN",n=nrow(y))
else{  
 ll<-list(coefficients=P[,1],logLik=logvero,AIC=AIC,BIC=BIC,iterations=iter,time=tempo,conv=conv,dist="MSL",class="MSMN",n=nrow(y))
 ll$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
}
 object.out<-ll
 class(object.out) <- "skewMLRM"
 object.out$y<-y.or
 object.out$X<-X.or
 object.out$"function"<-"estimate.MSL"
object.out
}
