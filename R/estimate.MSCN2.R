estimate.MSCN2 <-
function(y,X,max.iter=1000,prec=1e-4,est.var=TRUE,nu.fixed=0.1,gamma.fixed=0.5){
y.or<-y;X.or<-X
logscnnu<-function(nugama,res,A,Sigma){
nu=nugama[1]
gama=nugama[2]
fy=2*nu*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)+2*(1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*pnorm(A)
return(sum(log(fy)))}
logscngama<-function(gama,p,d,V,A){
fy=p/2*sum(V)*log(gama)-1/2*gama*sum(d*V)+sum(V*log(pnorm(sqrt(gama)*A)))
return(fy)}
mscn.logL <- function(theta,Y,X){
 n=nrow(Y)
 p=ncol(Y)
 q0=ncol(X[[1]])
 pth=length(theta)-2  # especifico
 beta0=as.matrix(theta[1:q0],q0,1)
 p1=p*(p+1)/2
 B= xpnd(theta[(q0+1):(q0+p1)]) # Sigma^{1/2}
 Sigma=B%*%B
 invSigma=solve2(Sigma)
 Binv=matrix.sqrt(invSigma) # Sigma^{-1/2} = B^{-1}
 lambda=as.matrix(theta[(q0+p1+1):pth],p,1)
 nu=as.numeric(theta[pth+1])
 gama=as.numeric(theta[pth+2])
 res<-matrix(0,n,p)
 for(i in 1:n){res[i,]<-Y[i,]-X[[i]]%*%beta0}
A=as.vector(t(lambda)%*%solve2(B)%*%t(res))
fy=2*nu*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)+2*(1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*pnorm(A)
return(sum(log(fy)))}
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
if(!nu.fixed & !gamma.fixed){
  aa=system.time({
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b0<-matrix(0,q,q)
  b1<-matrix(0,q,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,] }
  beta0<-solve2(b0)%*%b1
  res<-matrix(0,n,p)
  for(i in 1:n){res[i,]<-y[i,]-X[[i]]%*%beta0}
  S<-cov(res)
  invS<-solve2(S)
  B<-matrix.sqrt(S)
  lambda<-as.matrix(moments::skewness(res))
  nugama=c(0.5,0.5)
  theta0<-c(as.vector(beta0),vech(B),as.vector(lambda),nugama)
  log0= mscn.logL(theta0,y,X)
  Sigma=B%*%B
  invSigma=solve2(Sigma)
  delta=lambda/sqrt(1+sum(as.vector(lambda)^2))
  Delta= matrix.sqrt(solve2(Sigma))%*%delta
  a1=Sigma-Delta%*%t(Delta)
  bb=eigen(a1)$values
  { if (sum(bb>0)==p) Gama=a1
  else Gama=Sigma}
  A=as.vector(t(lambda)%*%solve2(B)%*%t(res))
  L1<-rbind(1e-3,1e-3)
  L2<-rbind(0.999,0.999)
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
    nu=nugama[1]
    gama=nugama[2] 
    fy=2*nu*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)+2*(1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*pnorm(A) 
    u=2/fy*(nu*gama*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)+ (1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*pnorm(A))
    eta1=2/fy*(nu*sqrt(gama)*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*dnorm(sqrt(gama)*A)+ (1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*dnorm(A))
    MT2=1/(1+as.numeric(t(Delta)%*%Gamainv%*%Delta))
    muT=MT2*as.vector(t(Delta)%*%Gamainv%*%t(res))
    ut = u*muT+sqrt(MT2)*eta1
    ut2=u*(muT^2)+MT2+muT*sqrt(MT2)*eta1
    # V estima
    V= 2*nu*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)/fy
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
    invSigma=solve2(Sigma)
    B=matrix.sqrt(Sigma)
    Binv=solve2(B)
    delta=Binv%*%Delta
    lambda=delta/sqrt(1-as.numeric(t(Delta)%*%invSigma%*%Delta))
    A=as.vector(t(lambda)%*%solve2(B)%*%t(res))    
    #nugama<-optim(nugama,logscnnu,gr = NULL,res,A,Sigma,method='L-BFGS-B',lower=L1,upper=L2,control=list(fnscale=-1))$par 
    # CM
    nu=mean(V)
    gamaf<-optim(gama,logscngama,gr=NULL,p,d,V,A,method="L-BFGS-B",lower=0.0001,upper=0.9999,control=list(fnscale=-1))
    gama<-as.numeric(gamaf$par)
    nugama=c(nu,gama)
    #
    theta=c(as.vector(beta0),vech(B),as.vector(lambda),nugama)
    #print(theta)
    logL=mscn.logL(theta,y,X)
    criterio=abs(logL-log0)
    bb=sum(eigen(Gama)$values>0)-p
    theta0=theta
    log0=logL
 }
})}
if(!nu.fixed & is.numeric(gamma.fixed)){
  aa=system.time({
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b0<-matrix(0,q,q)
  b1<-matrix(0,q,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,] }
  beta0<-solve2(b0)%*%b1
  res<-matrix(0,n,p)
  for(i in 1:n){res[i,]<-y[i,]-X[[i]]%*%beta0}
  S<-cov(res)
  invS<-solve2(S)
  B<-matrix.sqrt(S)
  lambda<-as.matrix(moments::skewness(res))
  nugama=c(0.5,gamma.fixed)
  theta0<-c(as.vector(beta0),vech(B),as.vector(lambda),nugama)
  log0= mscn.logL(theta0,y,X)
  Sigma=B%*%B
  invSigma=solve2(Sigma)
  delta=lambda/sqrt(1+sum(as.vector(lambda)^2))
  Delta= matrix.sqrt(solve2(Sigma))%*%delta
  a1=Sigma-Delta%*%t(Delta)
  bb=eigen(a1)$values
  { if (sum(bb>0)==p) Gama=a1
  else Gama=Sigma}
  A=as.vector(t(lambda)%*%solve2(B)%*%t(res))
  L1<-rbind(1e-3,1e-3)
  L2<-rbind(0.999,0.999)
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
    nu=nugama[1]
    gama=nugama[2] 
    fy=2*nu*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)+2*(1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*pnorm(A) 
    u=2/fy*(nu*gama*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)+ (1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*pnorm(A))
    eta1=2/fy*(nu*sqrt(gama)*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*dnorm(sqrt(gama)*A)+ (1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*dnorm(A))
    MT2=1/(1+as.numeric(t(Delta)%*%Gamainv%*%Delta))
    muT=MT2*as.vector(t(Delta)%*%Gamainv%*%t(res))
    ut = u*muT+sqrt(MT2)*eta1
    ut2=u*(muT^2)+MT2+muT*sqrt(MT2)*eta1
    # V estima
    V= 2*nu*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)/fy
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
    invSigma=solve2(Sigma)
    B=matrix.sqrt(Sigma)
    Binv=solve2(B)
    delta=Binv%*%Delta
    lambda=delta/sqrt(1-as.numeric(t(Delta)%*%invSigma%*%Delta))
    A=as.vector(t(lambda)%*%solve2(B)%*%t(res))    
    nu=mean(V)
    nugama=c(nu,gama)
    theta=c(as.vector(beta0),vech(B),as.vector(lambda),nugama)
    logL=mscn.logL(theta,y,X)
    criterio=abs(logL-log0)
    bb=sum(eigen(Gama)$values>0)-p
    theta0=theta
    log0=logL
 }
})}
if(is.numeric(nu.fixed) & !gamma.fixed){
  aa=system.time({
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  #array(1,c(p,1,n))
  q=ncol(X[[n]])
  b0<-matrix(0,q,q)
  b1<-matrix(0,q,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,] }
  beta0<-solve2(b0)%*%b1
  res<-matrix(0,n,p)
  for(i in 1:n){res[i,]<-y[i,]-X[[i]]%*%beta0}
  S<-cov(res)
  invS<-solve2(S)
  B<-matrix.sqrt(S)
  lambda<-as.matrix(moments::skewness(res))
  nugama=c(nu.fixed,0.5)
  theta0<-c(as.vector(beta0),vech(B),as.vector(lambda),nugama)
  log0= mscn.logL(theta0,y,X)
  Sigma=B%*%B
  invSigma=solve2(Sigma)
  delta=lambda/sqrt(1+sum(as.vector(lambda)^2))
  Delta= matrix.sqrt(solve2(Sigma))%*%delta
  a1=Sigma-Delta%*%t(Delta)
  bb=eigen(a1)$values
  { if (sum(bb>0)==p) Gama=a1
  else Gama=Sigma}
  A=as.vector(t(lambda)%*%solve2(B)%*%t(res))
  L1<-rbind(1e-3,1e-3)
  L2<-rbind(0.999,0.999)
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
    nu=nugama[1]
    gama=nugama[2] 
    fy=2*nu*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)+2*(1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*pnorm(A) 
    u=2/fy*(nu*gama*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)+ (1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*pnorm(A))
    eta1=2/fy*(nu*sqrt(gama)*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*dnorm(sqrt(gama)*A)+ (1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*dnorm(A))
    MT2=1/(1+as.numeric(t(Delta)%*%Gamainv%*%Delta))
    muT=MT2*as.vector(t(Delta)%*%Gamainv%*%t(res))
    ut = u*muT+sqrt(MT2)*eta1
    ut2=u*(muT^2)+MT2+muT*sqrt(MT2)*eta1
    # V estima
    V= 2*nu*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)/fy
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
    invSigma=solve2(Sigma)
    B=matrix.sqrt(Sigma)
    Binv=solve2(B)
    delta=Binv%*%Delta
    lambda=delta/sqrt(1-as.numeric(t(Delta)%*%invSigma%*%Delta))
    A=as.vector(t(lambda)%*%solve2(B)%*%t(res))    
    gamaf<-optim(gama,logscngama,gr=NULL,p,d,V,A,method="L-BFGS-B",lower=0.0001,upper=0.9999,control=list(fnscale=-1))
    gama<-as.numeric(gamaf$par)
    nugama=c(nu,gama)
    #
    theta=c(as.vector(beta0),vech(B),as.vector(lambda),nugama)
    #print(theta)
    logL=mscn.logL(theta,y,X)
    criterio=abs(logL-log0)
    bb=sum(eigen(Gama)$values>0)-p
    theta0=theta
    log0=logL
 }
})}
if(is.numeric(nu.fixed) & is.numeric(gamma.fixed)){
  aa=system.time({
  n=nrow(y)
  p=ncol(y)
  if(missing(X)) {X<-array(1,c(p,1,n))}
  q=ncol(X[[n]])
  b0<-matrix(0,q,q)
  b1<-matrix(0,q,1)
  for(i in 1:n){
    b0<-b0+t(X[[i]])%*%X[[i]]
    b1<-b1+t(X[[i]])%*%y[i,] }
  beta0<-solve2(b0)%*%b1
  res<-matrix(0,n,p)
  for(i in 1:n){res[i,]<-y[i,]-X[[i]]%*%beta0}
  S<-cov(res)
  invS<-solve2(S)
  B<-matrix.sqrt(S)
  lambda<-as.matrix(moments::skewness(res))
  nugama=c(nu.fixed,gamma.fixed)
  theta0<-c(as.vector(beta0),vech(B),as.vector(lambda),nugama)
  log0= mscn.logL(theta0,y,X)
  Sigma=B%*%B
  invSigma=solve2(Sigma)
  delta=lambda/sqrt(1+sum(as.vector(lambda)^2))
  Delta= matrix.sqrt(solve2(Sigma))%*%delta
  a1=Sigma-Delta%*%t(Delta)
  bb=eigen(a1)$values
  { if (sum(bb>0)==p) Gama=a1
  else Gama=Sigma}
  A=as.vector(t(lambda)%*%solve2(B)%*%t(res))
  L1<-rbind(1e-3,1e-3)
  L2<-rbind(0.999,0.999)
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
    nu=nugama[1]
    gama=nugama[2] 
    fy=2*nu*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)+2*(1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*pnorm(A) 
    u=2/fy*(nu*gama*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)+ (1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*pnorm(A))
    eta1=2/fy*(nu*sqrt(gama)*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*dnorm(sqrt(gama)*A)+ (1-nu)*mvtnorm::dmvnorm(res, sigma = Sigma)*dnorm(A))
    MT2=1/(1+as.numeric(t(Delta)%*%Gamainv%*%Delta))
    muT=MT2*as.vector(t(Delta)%*%Gamainv%*%t(res))
    ut = u*muT+sqrt(MT2)*eta1
    ut2=u*(muT^2)+MT2+muT*sqrt(MT2)*eta1
    # V estima
    V= 2*nu*mvtnorm::dmvnorm(res, sigma = Sigma/gama)*pnorm(sqrt(gama)*A)/fy
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
    invSigma=solve2(Sigma)
    B=matrix.sqrt(Sigma)
    Binv=solve2(B)
    delta=Binv%*%Delta
    lambda=delta/sqrt(1-as.numeric(t(Delta)%*%invSigma%*%Delta))
    A=as.vector(t(lambda)%*%solve2(B)%*%t(res))    
    theta=c(as.vector(beta0),vech(B),as.vector(lambda),nugama)
    logL=mscn.logL(theta,y,X)
    criterio=abs(logL-log0)
    bb=sum(eigen(Gama)$values>0)-p
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
 P<-matrix(theta[-(length(theta)-1:0)],ncol=1)
 nu<-theta[length(theta)-1]
 gamma<-theta[length(theta)]
 rownames(P)<-c(paste("beta",1:q,sep=""),paste("alpha",indices,sep=""),paste("lambda",1:p,sep=""))
 colnames(P)<-c("estimate") 
 conv.problem=1
 if(est.var)
 {
 MI.obs<-  FI.MSCN2(P,y,X,nu,gamma)
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
if(conv.problem==0) object.out<-list(coefficients=P[,1],se=P[,2],nu=nu,gamma=gamma,logLik=logL,AIC=AIC,BIC=BIC,iterations=cont,time=tempo,conv=conv,dist="MSCN2",class="MSMSN",n=nrow(y))
else
{
object.out<-list(coefficients=P[,1],nu=nu,gamma=gamma,logLik=logL,AIC=AIC,BIC=BIC,iterations=cont,time=tempo,conv=conv,dist="MSCN2",class="MSMSN",n=nrow(y))
object.out$warnings="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
}
 class(object.out) <- "skewMLRM"
 object.out$y<-y.or
 object.out$X<-X.or
 object.out$"function"<-"estimate.MSCN2"
 object.out
 }
