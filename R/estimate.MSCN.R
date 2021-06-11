estimate.MSCN <-
function(y,X,max.iter=1000,prec=1e-4,est.var=TRUE)
{
lmscr<-function(P,y,X){
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p*(p+1)/2)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):(length(P)-2)])
  v=as.numeric(P[length(P)-1])
  g=as.numeric(P[length(P)])
  e<-matrix(0,n,p)
  for(i in 1:n){e[i,]<-y[i,]-X[[i]]%*%b}
  aux=pmax(as.vector(t(h)%*%invB%*%t(e)),-37)
  l=sum(log(2*(v/det(1/sqrt(g)*B)*mvtnorm::dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*mvtnorm::dmvnorm(e%*%invB))*pnorm(aux)))
  return(l)}
lmscvgr<-function(PA,b,B,h,y,X){
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  v=PA[1]
  g=PA[2]
  invB=solve(B)
  e<-matrix(0,n,p)
  for(i in 1:n){e[i,]<-y[i,]-X[[i]]%*%as.matrix(b)}
  lvg=sum(log(v/det(1/sqrt(g)*B)*mvtnorm::dmvnorm(sqrt(g)*e%*%invB)+(1-v)/det(B)*mvtnorm::dmvnorm(e%*%invB)))
  return(lvg)}
escmscr<-function(P,y,X){
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  p0=p*(p+1)/2
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):(length(P)-2)])
  v=as.numeric(P[length(P)-1])
  g=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  W<-matrix(0,n,1)
  K<-matrix(0,n,1)
  u<-matrix(0,n,1)
  dldb1<-matrix(0,q,1)
  dldb2<-matrix(0,q,1)
  dlda<-matrix(0,p0,1)
  dldv<-0
  dldg<-0
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    K[i]<-sqrt(2*pi)*(v*g^(p/2)*dnorm(sqrt(g*d[i]))+(1-v)*dnorm(sqrt(d[i])))
    u[i]<-sqrt(2*pi)*(v*g^(p/2+1)*dnorm(sqrt(g*d[i]))+(1-v)*dnorm(sqrt(d[i])))/K[i]
    auxi<-pmax(as.vector(t(h)%*%invB%*%e[i,]),-37)
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
    dldv<-dldv+sqrt(2*pi)*(g^(p/2)*dnorm(sqrt(g*d[i]))-dnorm(sqrt(d[i])))/K[i]
    dldg<-dldg+sqrt(2*pi)*((p-g*d[i])*v*g^(p/2-1)*dnorm(sqrt(g*d[i])))/(2*K[i])}
  dldb=dldb1-dldb2
  dldh=invB%*%t(e)%*%W
  el=c(as.vector(dldb),as.vector(dlda),as.vector(dldh),dldv,dldg)
  return(el)
} 
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
  #t0=Sys.time()
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
  h<-as.matrix(moments::skewness(e))
  v<-0.5
  g<-0.5
  P_0<-c(as.vector(b),vech(B),as.vector(h),v,g)
  log0=lmscr(P_0,y,X)
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
    t<-e%*%invB%*%h+W #passo E
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
    g<-n*p*v/sum(v1*d) #passo M  
    P<-c(as.vector(b),vech(B),as.vector(h),v,g)
    logvero=lmscr(P,y,X)
    #crit<-sqrt(sum((P-P_0)^2))
    crit=abs(logvero-log0)
    P_0<-P
    log0=logvero
  }# o algoritmo EM
  })
  npar=length(P)
  AIC=-2*logvero+2*npar
  BIC=-2*logvero+log(nrow(y))*npar
  #t1=Sys.time()
  #escore=0#escmscr(P,y,X)
  #informacao=0#mifmscr(P,y,X)
 P<-matrix(P,ncol=1)
 colnames(P)<-c("estimate") 
 conv.problem=0
 if(est.var)
 {
 MI.obs<- FI.MSCN(P,y,X)
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
conv<-ifelse(iter<=max.iter & crit<=prec, 0, 1)
  tempo=as.numeric(aa[3])#as.numeric(t1-t0)
 aux=as.list(sapply(1:p,seq,by=1,to=p))
 indices=c()
 for(j in 1:p)
 {indices=c(indices,paste(j,aux[[j]],sep=""))}
 rownames(P)<-c(paste("beta",1:q,sep=""),paste("alpha",indices,sep=""),paste("lambda",1:p,sep=""),"nu","gamma")
  ll<-list(estimate=P,logLik=logvero,AIC=AIC,BIC=BIC,iterations=iter,time=tempo,conv=conv)
 if(conv.problem==1) ll$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
 object.out<-ll
object.out
}
