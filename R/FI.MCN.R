FI.MCN <-
function(P,y,X)
{
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
  v=as.numeric(P[length(P)-1])
  g=as.numeric(P[length(P)])
  d<-matrix(0,n,1)
  e<-matrix(0,n,p)
  K<-matrix(0,n,1)
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
  d2ldbdg<-matrix(0,q,1)
  d2ldadg<-matrix(0,p0,1)
  d2ldv2<-0
  d2ldvdg<-0
  d2ldg2<-0
  for(i in 1:n){
    e[i,]<-y[i,]-X[[i]]%*%b
    d[i]<-t(e[i,])%*%invB%*%invB%*%e[i,]
    K[i]<-sqrt(2*pi)*(v*g^(p/2)*dnorm(sqrt(g*d[i]))+(1-v)*dnorm(sqrt(d[i])))
    u[i]<-sqrt(2*pi)*(v*g^(p/2+1)*dnorm(sqrt(g*d[i]))+(1-v)*dnorm(sqrt(d[i])))/K[i]
    vu[i]<-sqrt(2*pi)*(v*(1-v)*(1-g)^2*g^(p/2)*dnorm(sqrt((g+1)*d[i])))/K[i]^2
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
    d2ldbdv<-d2ldbdv-sqrt(2*pi)*((1-g)*g^(p/2)*dnorm(sqrt((g+1)*d[i])))/K[i]^2*t(X[[i]])%*%invB%*%invB%*%e[i,]
    d2ldadv<-d2ldadv+sqrt(2*pi)*((1-g)*g^(p/2)*dnorm(sqrt((g+1)*d[i])))/(2*K[i]^2)*ddida
    d2ldbdg<-d2ldbdg-sqrt(2*pi)*v*((1-v)*(p-g*(2+p+d[i])+g^2*d[i])*g^(p/2-1)*dnorm(sqrt((g+1)*d[i]))-2*v*g^p*dnorm(sqrt(2*g*d[i])))/(2*K[i]^2)*t(X[[i]])%*%invB%*%invB%*%e[i,]
    d2ldadg<-d2ldadg+sqrt(2*pi)*v*((1-v)*(p-g*(2+p+d[i])+g^2*d[i])*g^(p/2-1)*dnorm(sqrt((g+1)*d[i]))-2*v*g^p*dnorm(sqrt(2*g*d[i])))/(4*K[i]^2)*ddida
    d2ldv2<-d2ldv2-(sqrt(2*pi)*(g^(p/2)*dnorm(sqrt(g*d[i]))-dnorm(sqrt(d[i])))/K[i])^2
    d2ldvdg<-d2ldvdg+sqrt(2*pi)*(p-g*d[i])*g^(p/2-1)*dnorm(sqrt((g+1)*d[i]))/(2*K[i]^2)  
    d2ldg2<-d2ldg2-sqrt(2*pi)*v*(2*p*v*g^(p-2)*dnorm(sqrt(2*g*d[i]))-(1-v)*((p-g*d[i])^2-2*p)*g^(p/2-2)*dnorm(sqrt((g+1)*d[i])))/(4*K[i]^2)
  }  
  d2ldb=cbind(d2ldb2,d2ldbda,d2ldbdh,d2ldbdv,d2ldbdg)
  d2lda=cbind(t(d2ldbda),d2lda2,t(d2ldhda),d2ldadv,d2ldadg)
  d2ldh=cbind(t(d2ldbdh),d2ldhda,d2ldh2,matrix(0,p,2))
  d2ldv=cbind(t(d2ldbdv),t(d2ldadv),matrix(0,1,p),as.matrix(d2ldv2),as.matrix(d2ldvdg))
  d2ldg=cbind(t(d2ldbdg),t(d2ldadg),matrix(0,1,p),as.matrix(d2ldvdg),as.matrix(d2ldg2))
  Hl=rbind(d2ldb,d2lda,d2ldh,d2ldv,d2ldg)
  return(Hl[-c(q+p0+1:p),-c(q+p0+1:p)])
}
