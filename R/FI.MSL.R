FI.MSL <-
function(P,y,X)
{
fauxs<-function(u,v,p,d,r,s){u^(v+p/2-r)*exp(-u*d/2)*log(u)^s}
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
  n=nrow(y)
  p=ncol(y)
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  b=as.matrix(P[1:q])
  B=xpnd(P[(q+1):(q+p0)])
  invB=solve(B)
  h=matrix(0, nrow=p)
  v=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  u<-matrix(0,n,1)
  vu<-matrix(0,n,1)
  W<-matrix(0,n,1)
  Wl<-matrix(0,n,1)
  d2ldb2<-matrix(0,q,q)
  d2ldbda<-matrix(0,q,p0)
  d2ldbdh<-matrix(0,q,p)
  d2lda2<-matrix(0,p0,p0)
  d2ldhda<-matrix(0,p,p0)
  d2ldh2<-matrix(0,p,p)
  d2ldbdv<-matrix(0,q,1)
  d2ldadv<-matrix(0,p0,1)
  d2ldv2<-matrix(0,1,1)
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    u[i]<-(2*v+p)/d[i]*pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i])
    vu[i]<-(2*v+p+2)*(2*v+p)/d[i]^2*pgamma(1,v+2+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i])-((2*v+p)/d[i])^2*(pgamma(1,v+1+p/2,scale=2/d[i])/pgamma(1,v+p/2,scale=2/d[i]))^2
    auxi<-pmax(t(h)%*%invB%*%e[i,],-37)
    W[i]<-dnorm(auxi)/pnorm(auxi)
    Wl[i]<--W[i]*(auxi+W[i])
    d2ldb2<-d2ldb2-as.numeric(u[i])*t(X[[i]])%*%invB%*%invB%*%X[[i]]+as.numeric(vu[i])*t(X[[i]])%*%invB%*%invB%*%e[i,]%*%t(e[i,])%*%invB%*%invB%*%X[[i]]+as.numeric(Wl[i])*t(X[[i]])%*%invB%*%h%*%t(h)%*%invB%*%X[[i]]
    d2ldbdh<-d2ldbdh-as.numeric(W[i])*t(X[[i]])%*%invB-as.numeric(Wl[i])*t(X[[i]])%*%invB%*%h%*%t(e[i,])%*%invB
    d2ldh2<-d2ldh2+as.numeric(Wl[i])*invB%*%e[i,]%*%t(e[i,])%*%invB
    ddida<-matrix(0,p0,1)
    dAida<-matrix(0,p0,1)
    d2Aidhda<-matrix(0,p,p0) 
    d2didbda<-matrix(0,q,p0)
    d2Aidbda<-matrix(0,q,p0)
    d2ldda2<-matrix(0,p0,p0)
    d2dida2<-matrix(0,p0,p0)
    d2Aida2<-matrix(0,p0,p0)
    for(j in 1:p0){
      ind=rep(0,p0)
      ind[j]=1
      Bj=xpnd(ind)
      d2Aidhda[,j]=-invB%*%Bj%*%invB%*%e[i,]
      dAida[j]=-t(h)%*%invB%*%Bj%*%invB%*%e[i,]
      d2Aidbda[,j]=t(X[[i]])%*%invB%*%Bj%*%invB%*%h
      ddida[j]=-t(e[i,])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
      d2didbda[,j]=2*t(X[[i]])%*%invB%*%(Bj%*%invB+invB%*%Bj)%*%invB%*%e[i,]
      for(k in 1:p0){
        ind=rep(0,p0)
        ind[k]=1
        Bk=xpnd(ind)
        d2ldda2[j,k]=-sum(diag(invB%*%Bk%*%invB%*%Bj))
        d2Aida2[j,k]=-t(h)%*%invB%*%(Bk%*%invB%*%Bj+Bj%*%invB%*%Bk)%*%invB%*%e[i,]
        d2dida2[j,k]=t(e[i,])%*%invB%*%(Bk%*%invB%*%Bj%*%invB+Bj%*%invB%*%Bk%*%invB+Bj%*%invB%*%invB%*%Bk+Bk%*%invB%*%invB%*%Bj+invB%*%Bk%*%invB%*%Bj+invB%*%Bj%*%invB%*%Bk)%*%invB%*%e[i,]
      }
    }
    d2ldhda<-d2ldhda+as.numeric(W[i])*d2Aidhda+as.numeric(Wl[i])*invB%*%e[i,]%*%t(dAida)
    d2ldbda<-d2ldbda-1/2*as.numeric(u[i])*d2didbda-1/2*as.numeric(vu[i])*t(X[[i]])%*%invB%*%invB%*%e[i,]%*%t(ddida)+as.numeric(W[i])*d2Aidbda-as.numeric(Wl[i])*t(X[[i]])%*%invB%*%h%*%t(dAida)
    d2lda2<-d2lda2-1/2*d2ldda2-1/2*as.numeric(u[i])*d2dida2+1/4*as.numeric(vu[i])*ddida%*%t(ddida)+as.numeric(W[i])*d2Aida2+as.numeric(Wl[i])*dAida%*%t(dAida)
    d2ldbdv<-d2ldbdv+(integrate(fauxs,lower=0,upper=1,v,p,d[i],0,1)$val/integrate(fauxs,lower=0,upper=1,v,p,d[i],1,0)$val-integrate(fauxs,lower=0,upper=1,v,p,d[i],1,1)$val*integrate(fauxs,lower=0,upper=1,v,p,d[i],0,0)$val/integrate(fauxs,lower=0,upper=1,v,p,d[i],1,0)$val^2)*t(X[[i]])%*%invB%*%invB%*%e[i,]
    d2ldadv<-d2ldadv-1/2*(integrate(fauxs,lower=0,upper=1,v,p,d[i],0,1)$val/integrate(fauxs,lower=0,upper=1,v,p,d[i],1,0)$val-integrate(fauxs,lower=0,upper=1,v,p,d[i],1,1)$val*integrate(fauxs,lower=0,upper=1,v,p,d[i],0,0)$val/integrate(fauxs,lower=0,upper=1,v,p,d[i],1,0)$val^2)*ddida
    d2ldv2<-d2ldv2+integrate(fauxs,lower=0,upper=1,v,p,d[i],1,2)$val/integrate(fauxs,lower=0,upper=1,v,p,d[i],1,0)$val-(integrate(fauxs,lower=0,upper=1,v,p,d[i],1,1)$val/integrate(fauxs,lower=0,upper=1,v,p,d[i],1,0)$val)^2-1/v^2
  }  
  d2ldb=cbind(d2ldb2,d2ldbda,d2ldbdh,d2ldbdv)
  d2lda=cbind(t(d2ldbda),d2lda2,t(d2ldhda),d2ldadv)
  d2ldh=cbind(t(d2ldbdh),d2ldhda,d2ldh2,matrix(0,p,1))
  d2ldv=cbind(t(d2ldbdv),t(d2ldadv),matrix(0,1,p),d2ldv2)
  Hl=rbind(d2ldb,d2lda,d2ldh,d2ldv)
  return(Hl[-c(q+p0+1:p),-c(q+p0+1:p)])
}
