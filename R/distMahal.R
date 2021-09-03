distMahal <-
function(object, alpha=0.95, ...)
{
qcn_alpha<-function(qalpha,alpha,p,nu,gama){
fy=nu*pchisq(qalpha*gama,p)+(1-nu)*pchisq(qalpha,p)
res=abs(fy-alpha)
return(res)}
qcn<-function(alpha,p,nu,gama){
upper1=qchisq(0.999,p)
q_alpha=optimize(qcn_alpha, c(0,upper1), alpha=alpha,p=p,nu=nu,gama=gama, lower = 0, upper = upper1)$minimum
return(as.numeric(q_alpha))
}
qsl_alpha<-function(qalpha,alpha,p,nu){
fy=pchisq(qalpha,p)-(2/qalpha)^nu*gamma(nu+p/2)/gamma(p/2)*pchisq(qalpha,p+2*nu)
res=abs(fy-alpha)
return(res)}
qsl<-function(alpha,p,nu){
upper1=qchisq(0.999,p+2*nu)
q_alpha=optimize(qsl_alpha, c(0,upper1), alpha=alpha,p=p,nu=nu, lower = 0, upper = upper1)$minimum
return(as.numeric(q_alpha))
}
dist=object$dist
y=object$y; X=object$X
  if (!any(dist == c("MN","MT","MSL","MCN","MSN","MSNC","MSTEC","MSTN","MSTT",
"MSSLEC","MSSL","MSSL2","MSCN","MSCN2","MSCEC")))
     stop("distribution is not recognized")
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
  theta=as.numeric(object$coefficients)
  n=nrow(y)
  q0=ncol(y)
  q1=q0*(q0+1)/2
  p=ncol(X[[n]])
  beta0<-theta[1:p]
  B<-xpnd(theta[(p+1):(p+q1)])
  Sigma=B%*%B
  Sigmainv=solve(Sigma)
  if(dist=="MT" | dist=="MSL") {nu=theta[p+q1+1]}
  if(dist=="MSTN" | dist=="MSSL") {nu=theta[p+q1+q0+1]}
  if(dist=="MCN") {nu=theta[p+q1+1];gamma=theta[p+q1+2]}
  if(dist=="MSCN") {nu=theta[p+q1+q0+1];gamma=theta[p+q1+q0+2]}
  if(dist=="MSTEC" | dist=="MSSLEC"| dist=="MSSL2"|dist=="MSTT") {nu=object$nu}
  if(dist=="MSCEC"|dist=="MSCN2") {nu=object$nu; gamma=object$gamma}
  Mahal=rep(0,n)
  for(i in 1:n){
      resi=matrix(y[i,]-X[[i]]%*%beta0,nrow=q0)
      Mahal[i]=t(resi)%*%Sigmainv%*%resi
      }
 pc=qchisq(alpha,q0)
 if(dist=="MT" | dist=="MSTN"|dist=="MSTT" | dist=="MSTEC" )  {  pc=qf(alpha,q0,nu)*q0 }
  if(dist=="MSL" | dist=="MSSL"|dist=="MSSL2" | dist=="MSSLEC" ) { pc=qsl(alpha,p,nu)
 }
  if(dist=="MCN" | dist=="MSCN"|dist=="MSCN2" | dist=="MSCEC" ) { pc=qcn(alpha,p,nu,gamma)
 }
 RVAL<-list(Mahal=Mahal)
 class(RVAL)<- "skewMLRM"
 RVAL$"function"<-"distMahal"
 RVAL$dist<-object$dist
 RVAL$class<-object$class
 RVAL$alpha<-alpha
 RVAL$cut<-pc
 RVAL$y<-object$y
 RVAL$X<-object$X
 return(RVAL)
}
