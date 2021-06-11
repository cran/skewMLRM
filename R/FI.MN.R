FI.MN <-
function(P,y,X){
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
 p=ncol(y)
 n=nrow(y)
 m=ncol(X[[1]])
 pth=length(P)
 beta=P[1:m]
 p1=p*(p+1)/2
 vecSigmam=P[(m+1):pth]
 B= xpnd(vecSigmam) # Sigma^{1/2}
 Sigma=B%*%B
 invSigma=solve(Sigma)
 Binv=matrix.sqrt(invSigma) # Sigma^{-1/2} = B^{-1}

 MI=matrix(0,pth,pth)
     
 for (i in 1:n){  
  resi=as.matrix(y[i,]-X[[i]]%*%beta)  # 
  ddimu=-2*invSigma%*%resi
  d2dimumu =2*t(X[[i]])%*%invSigma%*%X[[i]]
  #
  d2elltheta=matrix(0,pth,pth)
  ddialpha=matrix(0,p1,1)
  d2dialphamu=matrix(0,m,p1)
  d2dialpha=matrix(0,p1,p1)
  d2logBalpha=matrix(0,p1,p1)
  #
  d2logtheta=matrix(0,pth,pth)
  d2ditheta=matrix(0,pth,pth)
 
   
  for (k in 1:p1){
        ind=rep(0,p1)
        ind[k]=1
        Bk=xpnd(ind)
        #ddialpha[k]=-t(resi)%*%Binv%*%(Bk%*%Binv+Binv%*%Bk)%*%Binv%*%resi
        d2dialphamu[,k]=2*t(X[[i]])%*%Binv%*%(Bk%*%Binv+Binv%*%Bk)%*%Binv%*%resi
        for (s in 1:p1){
            ind=rep(0,p1)
            ind[s]=1
            Bs=xpnd(ind)
            d2logBalpha[k,s]=-sum(diag(Binv%*%Bs%*%Binv%*%Bk))
            d2dialpha[k,s]=t(resi)%*%Binv%*%(Bs%*%Binv%*%Bk%*%Binv+Bk%*%Binv%*%Bs%*%Binv+Bk%*%invSigma%*%Bs+Bs%*%invSigma%*%Bk+Binv%*%Bs%*%Binv%*%Bk+Binv%*%Bk%*%Binv%*%Bs)%*%Binv%*%resi
        }
    }
    #
    d2logtheta[(m+1):pth,(m+1):pth]= d2logBalpha
    #
    d2ditheta[1:m,1:m]= d2dimumu
    d2ditheta[1:m,(m+1):pth]= d2dialphamu
    d2ditheta[(m+1):pth,(m+1):pth]= d2dialpha
    d2ditheta[(m+1):pth,1:m]= t(d2dialphamu)      
    #
    d2elltheta=-1/2*d2logtheta-1/2*d2ditheta
    MI=MI+d2elltheta    
 }
 #print(diag(solve(-MI)))
 return(MI)
 }
