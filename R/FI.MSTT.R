FI.MSTT <-
function(P,y,X,nu)
{
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
  n=nrow(y)
  p=ncol(y)
  q=ncol(X[[n]])
  p0=p*(p+1)/2
  beta=as.matrix(P[1:q])
  Dr=xpnd(P[(q+1):(q+p0)])
  Sigma=Dr%*%Dr
  Dr.inv=solve(Dr)
  d.sig <- det(Sigma) 
  lambda=as.matrix(P[(q+p0+1):length(P)])
  nu=as.numeric(nu)
  I.Phi <- function(w=0,Ai=NULL,di,nu=0) {as.numeric((( 2^w*nu^(nu/2)*gamma(w + nu/2))/(gamma(nu/2)*(nu + di)^(nu/2 + w)))*pt( ((Ai)/(di + nu)^(0.5))*sqrt(2*w + nu), 2*w + nu))}
  I.phi <- function(w=0,Ai=NULL,di,nu=0) {as.numeric(((2^w*nu^(nu/2))/(sqrt(2*pi)*gamma(nu/2)))*(1/(di + Ai^2 + nu))^((nu + 2*w)/2)*gamma((nu + 2*w)/2))}
  deriv.der <- function(A,B,C) det(A)*sum(B * t(C))
npar=q+p*(p+1)/2 + p    # beta + Sigma + lambda
soma<-matrix(0,npar,npar)
for (i in 1:n){
yi <- matrix(y[i,], p, 1)
Xi<- X[[i]] 
mui<-matrix(Xi%*%beta,p,1)
Ai <- as.numeric(t(lambda)%*%Dr.inv%*%(yi - mui))
di <- mahalanobis(as.vector(yi), as.vector(mui), Sigma)
dAir.dbeta <- -t(Xi)%*%Dr.inv%*%lambda
dAir.dlambda <- Dr.inv%*%(yi - mui)
dir.dbeta <- -2*t(Xi)%*%(Dr.inv%*%Dr.inv)%*%(yi - mui)
ddet.ds<-c()
dir.ds<-c()
dAir.ds<-c()
                    
l <- m <- 1
for(k in 1:((p+1)*p/2)) {
Vis <- FALSE
D <- matrix(rep(0,p*p),p,p)
D[l,m] <- D[m,l] <- 1
ddet.ds[k] <- -(1/det(Dr)^2)*deriv.der(Dr,Dr.inv,D)
dir.ds[k] <- - t(yi - mui)%*%Dr.inv%*%(D%*%Dr.inv + Dr.inv%*%D)%*%Dr.inv%*%(yi - mui)
dAir.ds[k] <- - t(lambda)%*%Dr.inv%*%D%*%Dr.inv%*%(yi - mui)

if(((l*m - p*floor((l*m)/p)) == 0) && (l != m)) {
l <- l+1
m <- l
Vis <- TRUE
                                                } 
           if(!Vis) m <- m+1
                          } 
ki=I.Phi(w=p/2,Ai=Ai,di=di,nu=nu)
derkibeta=I.phi(w=(p+1)/2,Ai=Ai,di=di,nu=nu)*dAir.dbeta-0.5*I.Phi(w=(p+2)/2,Ai=Ai,di=di,nu=nu)*dir.dbeta
Sbeta=(1/ki)*derkibeta
derkisigma=I.phi(w=(p+1)/2,Ai=Ai,di=di,nu=nu)*dAir.ds-0.5*I.Phi(w=(p+2)/2,Ai=Ai,di=di,nu=nu)*dir.ds
Ssigma=-0.5*ddet.ds+(1/ki)*derkisigma
derkilambda=I.phi(w=(p+1)/2,Ai=Ai,di=di,nu=nu)*dAir.dlambda
Slambda=(1/ki)*derkilambda     
S=matrix(c(Sbeta,Ssigma,Slambda),npar,1)
soma <- soma + S%*%t(S)
     }       
return(-soma)
}
