estimate.MN <-
function(y,X,max.iter=1000,prec=1e-4,est.var=TRUE)
{
lmnr<-function(theta,y,X){
  n=nrow(y)
  q0=ncol(y)
  if(missing(X)) {X<-array(1,c(q0,1,n))}
  p=ncol(X[[n]])
  beta0=as.matrix(theta[1:p])
  B=xpnd(theta[(p+1):(p+q0*(q0+1)/2)])
  Sigma=B%*%B
  l=0
    for(i in 1:n){
    mui=as.vector(X[[i]]%*%beta0)
    l=l+log(mvtnorm::dmvnorm(as.vector(y[i,]),mui,Sigma))}
  return(l)}
X.or<-X
y.or<-y
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
  t0=Sys.time()
 aa=system.time({
  n=nrow(y)
  q0=ncol(y)
  p=ncol(X[[n]])
  m=nrow(X[[n]])
  b0<-matrix(0,p,p)
  b1<-matrix(0,p,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,] 
  }
  beta0<-solve2(b0)%*%b1
  e<-matrix(0,n,q0)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%beta0
  }
  Sigma<-cov(e)
  invS<-solve2(Sigma)
  B<-matrix.sqrt(Sigma)
  invB<-solve2(B)
  P_0<-c(as.vector(beta0),vech(B))
  log0=lmnr(P_0,y,X)
  crit<-1
  iter<-0
  while((crit>=prec)&&(iter<=max.iter))
  { 
    iter<-iter+1 
    b0<-matrix(0,p,p)
    b1<-matrix(0,p,1)
    for(i in 1:n){
      b0<-b0+t(X[[i]])%*%invS%*%X[[i]]
      b1<-b1+t(X[[i]])%*%invS%*%y[i,]
    }
    beta0<-solve2(b0)%*%b1
    e<-matrix(0,n,m)
    for(i in 1:n){
      e[i,]<-y[i,]-X[[i]]%*%beta0}
    Sigma<-1/n*t(e)%*%e
    invS<-solve2(Sigma)
    B<-matrix.sqrt(Sigma)
    P<-c(as.vector(beta0),vech(B))
    log1=lmnr(P,y,X)
    crit=abs(log1-log0)
    P_0<-P
    log0=log1
  }
 })
  P=matrix(P,ncol=1)
  npar=length(P)
conv.problem=1
 if(est.var)
 {
 MI.obs<-  FI.MN(P,y,X)
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
conv<-ifelse(iter<=max.iter & crit<=prec, 0, 1)
   tempo=as.numeric(aa[3])
  logvero<-lmnr(P,y,X)  
  AIC=-2*logvero+2*npar
  BIC=-2*logvero+log(nrow(y))*npar
 aux=as.list(sapply(1:m,seq,by=1,to=m))
 tempo=as.numeric(aa[3])
 indices=c()
 for(j in 1:m)
 {indices=c(indices,paste(j,aux[[j]],sep=""))}
 rownames(P)<-c(paste("beta",1:p,sep=""),paste("alpha",indices,sep=""))
if(conv.problem==0)  ll<-list(coefficients=P[,1],se=P[,2],logLik=logvero,AIC=AIC,BIC=BIC,iterations=iter,time=tempo,conv=conv,dist="MN",class="MSMN",n=nrow(y))
else{
ll<-list(coefficients=P[,1],logLik=logvero,AIC=AIC,BIC=BIC,iterations=iter,time=tempo,conv=conv,dist="MN",class="MSMN",n=nrow(y));
ll$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"}
 object.out<-ll
 class(object.out) <- "skewMLRM"
 object.out$y<-y.or
 object.out$X<-X.or
 object.out$"function"<-"estimate.MN"
 object.out
}
