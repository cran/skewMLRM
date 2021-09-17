FI.MSNC <-
function(P,y,X)
{
          if (missing(X)) {
            X <- array(1, c(q0, 1, n))
        }
    y <- as.matrix(y)
    if (!is.matrix(y)) 
        stop("y must have at least one element")
    if (is.null(X)) {
        X <- array(c(diag(ncol(y))), c(ncol(y), ncol(y), nrow(y)))
    }
    if (is.array(X) == FALSE & is.list(X) == FALSE) 
        stop("X must be an array or a list")
    if (is.array(X)) {
        Xs <- list()
        for (i in 1:nrow(y)) {
            Xs[[i]] <- matrix(t(X[, , i]), nrow = ncol(y))
        }
        X <- Xs
    }
    if (ncol(y) != nrow(X[[1]])) 
        stop("y does not have the same number of columns than X")
    if (nrow(y) != length(X)) 
        stop("y does not have the same number of observations than X")
 # theta=list(beta=beta.new,Sigma=Sigma.new,eta=eta.new,dif=dif,iter=i)
  p=ncol(y)
  n=nrow(y)
 m=ncol(X[[1]])
 beta=P[1:m]
  q0=length(beta)
  Sigma<-matrix(0,ncol=p, nrow=p)
  Sigma[lower.tri(Sigma, diag = TRUE)]<-P[1:(p*(p+1)/2)+m]
  Sigma[upper.tri(Sigma, diag = TRUE)]<-P[1:(p*(p+1)/2)+m]
  Sinv=solve2(Sigma) # Sigma^{-1}
  C=matrix.sqrt(Sigma)
  Cinv=solve2(C)# Sigma^{-1/2}
  eta=P[m+(p*(p+1)/2)+1:p]
  lambda=C%*%eta
  q1=p*(p+1)/2
  pth=q0+p+q1

  dlogS=matrix(0,pth,1)
        # dlog Sigma
      for (r in 1:q1){
        indr=rep(0,q1)
        indr[r]=1
        Cr=xpnd(indr)
        auxi=Sinv%*%(Cr%*%C+C%*%Cr)
        dlogS[q0+r]=sum(diag(auxi))
        }

  dAi=matrix(0,pth,1)
  ddi=matrix(0,pth,1)
  MI=matrix(0,pth,pth)
  for (i in 1:n){
      resi=y[i,]-X[[i]]%*%beta  # p x 1
      ai= as.numeric(t(eta)%*%resi)
      di= as.numeric(t(resi)%*%Sinv%*%resi)
      Iphi.i<-exp(-0.5*di)*dcauchy(ai)
      Ki<-exp(-0.5*di)*pcauchy(ai)
      
      dAi[1:q0]=t(X[[i]])%*%Cinv%*%lambda
      dAi[(q0+q1+1):pth]=Cinv%*%resi
      ddi[1:q0]=-2*t(X[[i]])%*%Cinv%*%Cinv%*%resi
      for (r in 1:q1){
           indr=rep(0,q1)
           indr[r]=1
           Cr=xpnd(indr)
           ddi[q0+r]=-t(resi)%*%Cinv%*%(Cr%*%C+C%*%Cr)%*%Cinv%*%resi
           dAi[q0+r]=-t(lambda)%*%Cinv%*%Cr%*%Cinv%*%resi
           }
      scoretheta=-0.5*dlogS +(Iphi.i*dAi-0.5*Ki*ddi)/Ki
      MI=MI+scoretheta%*%t(scoretheta)
    } 
 return(-MI)
}
