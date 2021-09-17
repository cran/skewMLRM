choose.MSMN <-
function(y, X=NULL, max.iter=1000, prec=1e-4, est.var=TRUE, criteria="AIC",cluster=FALSE)
{

se.est<-function(P,y,X,dist="MN", nu=3, gamma=0.5)
{
#mf <- match.call(expand.dots = FALSE)
#if("nu" %in% names(mf$...)) nu=mf$...$nu
#if("gamma" %in% names(mf$...)) gamma=mf$...$gamma
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
  FI.dist <- get(paste("FI.", dist, sep = ""), mode = "function")
if(dist=="MSTT" | dist=="MSSL2" | dist=="MSTEC" | dist=="MSSLEC") {
#if(!exists("nu")) nu=3 
if(!is.numeric(nu)) stop("nu must be a number greater than 1") 
if(as.numeric(nu)<1) stop("nu must be a number greater than 1")}
if(dist=="MSCN2" | dist=="MSCEC") {
#if(!exists("nu")) nu=0.1 
#if(!exists("gamma")) gamma=0.5 
if(!is.numeric(nu)) stop("nu must be a number between 0 and 1") 
if((as.numeric(nu)<0 | as.numeric(nu)>1)) stop("nu must be a number between 0 and 1") 
if(!is.numeric(gamma)) stop("gamma must be a number between 0 and 1")
if(as.numeric(gamma)<0 | as.numeric(gamma)>1) stop("gamma must be a number between 0 and 1")}
P<-matrix(P,ncol=1)
if (!any(dist == c("MSTEC","MSSLEC","MSCEC","MSTT","MSSL2","MSCN2"))) MI.obs<-FI.dist(P=P, y=y, X=X)
if (dist == "MSTEC" | dist=="MSSLEC" | dist=="MSTT" | dist=="MSSL2") MI.obs<-FI.dist(P=P, y=y, X=X, nu=nu)
if (dist == "MSCEC" | dist=="MSCN2") MI.obs<-FI.dist(P=P, y=y, X=X, nu=nu, gamma=gamma)
 test=try(solve2(MI.obs),silent=TRUE)
 se=c()
 if(is.numeric(test) & max(diag(test))<0) 
 {
 se=sqrt(-diag(test))
 }
 else  stop("Standard errors can't be estimated: Numerical problems with the inversion of the information matrix")
 se
}
  if(cluster!=TRUE) cluster=FALSE
  if(all(criteria!=c("AIC","BIC")))
        stop("criteria must be AIC or BIC")
  if (cluster==FALSE){      
  fit.mn=try(estimate.MN(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mtn=try(estimate.MT(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.msl=try(estimate.MSL(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mcn=try(estimate.MCN(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  }
  else {
   cl0 <- parallel::detectCores()-1
   cl <- parallel::makeCluster(cl0)
  doParallel::registerDoParallel(cl)
      kname=c("MN","MT","MSL","MCN");j=1
  fit4 <- foreach(j = 1:4) %dopar%{
     #source('SMNfunctions.r')
     require(skewMLRM)
     XX=kname[j]
     estimate.dist.XX <- get(paste("estimate.", XX, sep = ""), mode = "function")
     fit4= try(estimate.dist.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)  
   }
   fit.mn=fit4[[1]]
   fit.mtn=fit4[[2]]
   fit.msl=fit4[[3]]
   fit.mcn=fit4[[4]]
  parallel::stopCluster(cl)
   }
  if(grepl("Error",fit.mn)[1] & grepl("Error",fit.mtn)[1] & 
     grepl("Error",fit.msl)[1] & grepl("Error",fit.mcn)[1])
  {stop("estimation problem in all the models in the SMN class")}
  index<-c(); maxi=c()
  if(criteria=="AIC")
  {
  if(grepl("Error",fit.mtn)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mtn$AIC)}
  if(grepl("Error",fit.msl)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.msl$AIC)}
  if(grepl("Error",fit.mcn)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mcn$AIC)}
  if(grepl("Error",fit.mn)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.mn$AIC)}
  }
  if(criteria=="BIC")
  {
  if(grepl("Error",fit.mtn)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mtn$BIC)}
  if(grepl("Error",fit.msl)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.msl$BIC)}
  if(grepl("Error",fit.mcn)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mcn$BIC)}
  if(grepl("Error",fit.mn)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.mn$BIC)}
  }
  fitted<-c("MT","MSL","MCN","MN")[index]
  index<-index[which.min(maxi)]
  fit<-switch(index, "1"=fit.mtn, "2"=fit.msl, "3"=fit.mcn, "4"=fit.mn)
  selected<-c("MT","MSL","MCN","MN")[index]
  se<-"Error"
  if(est.var)
  {
   se<-try(se.est(fit$coefficients,y,X,dist=selected),silent=TRUE)
  } 
  if(!grepl("Error",se[1])){
  names(se)<-names(fit$coefficients)
  RVAL<- list(coefficients=fit$coefficients, se=se, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMN", selected.model=selected)  
  }
  else{
  RVAL<- list(coefficients=fit$coefficients, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMN", selected.model=selected)
  RVAL$warning="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
}
 class(RVAL)<- "skewMLRM"
 RVAL$y<-y; RVAL$X<-X
 RVAL$dist<-RVAL$selected.model
 RVAL$choose.crit<-criteria
 RVAL$"function"<-"choose.MSMN"
 RVAL
}
