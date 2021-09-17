estimate.MSTN <-
function(y,X,max.iter=1000,prec=1e-4,est.var=TRUE,nu.min=2.0001)
{
y.or<-y; X.or<-X
lmstr<-function(P,y,X){
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p*(p+1)/2)])
  invB=solve2(B)
  h=as.matrix(P[(q+p*(p+1)/2+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]}
  aux=pmax(as.vector(t(h)%*%invB%*%t(e)),-37)
  l=sum(log(2*gamma((v+p)/2)/(gamma(v/2)*(pi*v)^(p/2)*det(B))*(1+d/v)^(-(v+p)/2)*pnorm(aux)))
  return(l)}
lmstvr<-function(v,b,B,h,y,X){
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  invB=solve2(B)
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%as.matrix(b)
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]}
  lv=sum(log(2*gamma((v+p)/2)/(gamma(v/2)*(pi*v)^(p/2))*(1+d/v)^(-(v+p)/2)))
  return(lv)} 
escmstr<-function(P,y,X){
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b=as.matrix(P[1:q])
  p0=p*(p+1)/2
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve2(B)
  h=as.matrix(P[(q+p0+1):(length(P)-1)])
  v=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  W<-matrix(0,n,1)
  u<-matrix(0,n,1)
  dldb1<-matrix(0,q,1)
  dldb2<-matrix(0,q,1)
  dlda<-matrix(0,p0,1)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(v+p)/(v+d[i])
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
      ddida[j]<--t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
    }
    dlda<-dlda-1/2*dldda-1/2*as.numeric(u[i])*ddida+as.numeric(W[i])*dAida
  }
  dldb=dldb1-dldb2
  dldh=invB%*%t(e)%*%W
  dldv=n/2*(digamma((v+p)/2)-log((v+p)/2)-digamma(v/2)+log(v/2)+1)-1/2*sum(u-log(u))
  el=c(as.vector(dldb),as.vector(dlda),as.vector(dldh),dldv)
  return(el)}
 y<-as.matrix(y)
if(!is.matrix(y)) stop("y must have at least one element")
if(!is.numeric(nu.min) | nu.min<=0) stop("nu.min should be a positive number")
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
  b<-solve2(b0)%*%b1
  e<-matrix(0,n,p)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
  }
  S<-cov(e)
  invS<-solve2(S)
  B<-matrix.sqrt(S)
  invB<-solve2(B)
  h<-as.matrix(moments::skewness(e))
  v<-5
  P_0<-c(as.vector(b),vech(B),as.vector(h),v)
  log0=lmstr(P_0,y,X)
  crit<-1
  iter<-0
  while((crit>=prec)&&(iter<=max.iter))
  { 
    iter<-iter+1 
    d<-matrix(0,n,1)
    e<-matrix(0,n,p)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%b
      d[i]<-t(e[i,])%*%invS%*%e[i,]
    }
    u<-(v+p)/(v+d)
    #lu<-digamma((v+p)/2)-log((v+d)/2)
    aux<-pmax(as.vector(e%*%invB%*%h),-37)
    W<-as.matrix(dnorm(aux)/pnorm(aux))
    t<-e%*%invB%*%h+W #passo E
    S<-1/n*(t(e)%*%diag(as.vector(u),n)%*%e)
    invS<-solve2(S)
    B<-matrix.sqrt(S)
    invB<-solve2(B)
    h<-B%*%solve2(t(e)%*%e)%*%t(e)%*%t
    b0<-matrix(0,q,q)
    b1<-matrix(0,q,1)
    for(i in 1:n){
      b0<-b0+t(X[[i]])%*%invB%*%(as.numeric(u[i])*diag(p)+h%*%t(h))%*%invB%*%X[[i]]
      b1<-b1+t(X[[i]])%*%(as.numeric(u[i])*(invB%*%invB)%*%y[i,]-as.numeric(t[i])*invB%*%h+invB%*%h%*%t(h)%*%invB%*%y[i,])
    }
    b<-solve2(b0)%*%b1
    V<-optim(v,lmstvr,gr=NULL,b,B,h,y,X,method="L-BFGS-B",lower=nu.min,upper=100,control=list(fnscale=-1,maxit=20))
    v<-as.numeric(V$par) #passo M
    P<-c(as.vector(b),vech(B),as.vector(h),v)
    logvero=lmstr(P,y,X)
    #crit<-sqrt(sum((P-P_0)^2))
    crit=abs(logvero-log0)
    P_0<-P
    log0=logvero
  }
  })
  npar=length(P)
  AIC=-2*logvero+2*npar
  BIC=-2*logvero+log(nrow(y))*npar
 P<-matrix(P,ncol=1)
 colnames(P)<-c("estimate") 
 conv.problem=1
 if(est.var){
 MI.obs<-  FI.MSTN(P,y,X)
 test=try(solve2(MI.obs),silent=TRUE)
 se=c()
 if(is.numeric(test) & max(diag(test))<0){
 conv.problem=0
 se=sqrt(-diag(test))
 P<-cbind(P,se)
 colnames(P)<-c("estimate","s.e.")}
}
conv<-ifelse(iter<=max.iter & crit<=prec, 0, 1)
 tempo=as.numeric(aa[3])#as.numeric(t1-t0)
 aux=as.list(sapply(1:p,seq,by=1,to=p))
 indices=c()
 for(j in 1:p)
 {indices=c(indices,paste(j,aux[[j]],sep=""))}
 rownames(P)<-c(paste("beta",1:q,sep=""),paste("alpha",indices,sep=""),paste("lambda",1:p,sep=""),"nu")
if(conv.problem==0) ll<-list(coefficients=P[,1],se=P[,2],logLik=logvero,AIC=AIC,BIC=BIC,iterations=iter,time=tempo,conv=conv,dist="MSTN",class="MSSMN",n=nrow(y))
else{
ll<-list(coefficients=P[,1],logLik=logvero,AIC=AIC,BIC=BIC,iterations=iter,time=tempo,conv=conv,dist="MSTN",class="MSSMN",n=nrow(y))
ll$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
}
 object.out<-ll
 class(object.out) <- "skewMLRM"
 object.out$y<-y.or
 object.out$X<-X.or
 object.out$"function"<-"estimate.MSTN"
 object.out
}
