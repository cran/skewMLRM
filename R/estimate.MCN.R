estimate.MCN <-
function(y,X,max.iter=1000,prec=1e-4,est.var=TRUE)
{
  y.or<-y;X.or<-X
lmcr<-function(P,y,X){
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p*(p+1)/2)])
  invB=solve(B)
  h=matrix(0, nrow=p)
  v=as.numeric(P[length(P)-1])
  g=as.numeric(P[length(P)])
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b}
  aux=pmax(as.vector(t(h)%*%invB%*%t(e)),-37)
  l=sum(log(2*(v/det(1/sqrt(g)*B)*mvtnorm::dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*mvtnorm::dmvnorm(e%*%invB))*pnorm(aux)))
  return(l)} 
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
  B<-matrix.sqrt(S)
  invB<-solve(B)
  h<-matrix(0, nrow=p)
  v<-0.5
  g<-0.5
  P_0<-c(as.vector(b),vech(B),as.vector(h),v,g)
  log0=lmcr(P_0,y,X)
  crit<-1
  iter<-0
  while((crit>=prec)&&(iter<=max.iter))
  { 
    iter<-iter+1 
    v1<-v/det(1/sqrt(g)*B)*mvtnorm::dmvnorm(sqrt(g)*e%*%invB)/(v/det(1/sqrt(g)*B)*mvtnorm::dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*mvtnorm::dmvnorm(e%*%invB))
    v2<-(1-v)/det(B)*mvtnorm::dmvnorm(e%*%invB)/(v/det(1/sqrt(g)*B)*mvtnorm::dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*mvtnorm::dmvnorm(e%*%invB))
    u<-g*v1+v2
    aux<-pmax(as.vector(e%*%invB%*%h),-37)
    W<-as.matrix(dnorm(aux)/pnorm(aux))
    t<-e%*%invB%*%h+W 
    S<-1/n*(t(e)%*%diag(as.vector(u))%*%e)
    invS<-solve(S)
    B<-matrix.sqrt(S)
    invB<-solve(B)
    h<-B%*%solve(t(e)%*%e)%*%t(e)%*%t
    b0<-matrix(0,q,q)
    b1<-matrix(0,q,1)
    for(i in 1:n){
      b0<-b0+t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
      b1<-b1+t(X[[i]])%*%(as.numeric(u[i])*invS%*%y[i,]-as.numeric(t[i])*invB%*%h+invB%*%h%*%t(h)%*%invB%*%y[i,])
    }
    C<-solve(b0)
    b<-C%*%b1
    d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b
      d[i]<-t(e[i,])%*%invS%*%e[i,]
    }
    v<-sum(v1)/n
    g<-n*p*v/sum(v1*d)   
    P<-c(as.vector(b),vech(B),as.vector(h),v,g)
    logvero=lmcr(P,y,X)
    crit=abs(logvero-log0)
    P_0<-P
    log0=logvero
  }
  })
  npar=length(P)
  AIC=-2*logvero+2*npar
  BIC=-2*logvero+log(nrow(y))*npar
 P<-matrix(P[-c(q+p*(p+1)/2+1:p)],ncol=1)
 colnames(P)<-c("estimate") 
 conv.problem=1
 if(est.var)
 {
 MI.obs<- FI.MCN(P,y,X)
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
conv<-ifelse(iter<=max.iter & crit<=prec, 0, 1)
  tempo=as.numeric(aa[3])
 aux=as.list(sapply(1:p,seq,by=1,to=p))
 indices=c()
 for(j in 1:p)
 {indices=c(indices,paste(j,aux[[j]],sep=""))}
 rownames(P)<-c(paste("beta",1:q,sep=""),paste("alpha",indices,sep=""),"nu","gamma")
if(conv.problem==0) ll<-list(coefficients=P[,1],se=P[,2],logLik=logvero,AIC=AIC,BIC=BIC,iterations=iter,time=tempo,conv=conv,dist="MCN",class="MSMN",n=nrow(y))
else{
ll<-list(coefficients=P[,1],logLik=logvero,AIC=AIC,BIC=BIC,iterations=iter,time=tempo,conv=conv,dist="MCN",class="MSMN",n=nrow(y))
ll$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
}
 object.out<-ll
 class(object.out) <- "skewMLRM"
 object.out$y<-y.or
 object.out$X<-X.or
 object.out$"function"<-"estimate.MCN"
object.out
}
