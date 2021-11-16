choose.MSMSN <-
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
  fit.msn=try(estimate.MSN(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mstt=try(estimate.MSTT(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE, nu.fixed = FALSE),silent=TRUE)
  fit.mssl2=try(estimate.MSSL2(y,X,max.iter=max.iter,prec=prec,est.var=FALSE, nu.fixed = FALSE),silent=TRUE)
  fit.mscn2=try(estimate.MSCN2(y,X,max.iter=max.iter,prec=prec,est.var=FALSE, nu.fixed = FALSE, gamma.fixed=FALSE),silent=TRUE)
  }
  else {
   cl0 <- parallel::detectCores()-1
   cl <- parallel::makeCluster(cl0)
  doParallel::registerDoParallel(cl);j=1
      kname=c("MSN","MSTT","MSSL2","MSCN2")
  fit4 <- foreach(j = 1:4) %dopar%{
     #source('SMSNfunctions.r')
     require(skewMLRM)
     XX=kname[j]
     estimate.dist.XX <- get(paste("estimate.", XX, sep = ""), mode = "function")
     if(XX=="MSN") fit4= try(estimate.dist.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)  
     if(XX=="MSTT" | XX=="MSSL2") fit4= try(estimate.dist.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE,nu.fixed=FALSE),silent=TRUE)  
     if(XX=="MSCN2") fit4= try(estimate.dist.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE,nu.fixed=FALSE,gamma.fixed=FALSE),silent=TRUE)  
   }
   fit.msn=fit4[[1]]
   fit.mstt=fit4[[2]]
   fit.mssl2=fit4[[3]]
   fit.mscn2=fit4[[4]]
  parallel::stopCluster(cl)
   }
  if(grepl("Error",fit.mstt)[1] & grepl("Error",fit.mssl2)[1] & 
grepl("Error",fit.mscn2)[1] & grepl("Error",fit.msn)[1])
  {stop("estimation problem in all the models in the SMSN class")}
  index<-c(); maxi=c()
  if(criteria=="AIC")
  {
  if(grepl("Error",fit.mstt)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mstt$AIC)}
  if(grepl("Error",fit.mssl2)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.mssl2$AIC)}
  if(grepl("Error",fit.mscn2)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mscn2$AIC)}
  if(grepl("Error",fit.msn)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msn$AIC)}
  }
  if(criteria=="BIC")
  {
  if(grepl("Error",fit.mstt)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mstt$BIC)}
  if(grepl("Error",fit.mssl2)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.mssl2$BIC)}
  if(grepl("Error",fit.mscn2)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mscn2$BIC)}
  if(grepl("Error",fit.msn)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msn$BIC)}
  }
  fitted<-c("MSTT","MSSL2","MSCN2","MSN")[index]
  index<-index[which.min(maxi)]
  fit<-switch(index, "1"=fit.mstt, "2"=fit.mssl2, "3"=fit.mscn2, "4"=fit.msn)
  selected<-c("MSTT","MSSL2","MSCN2","MSN")[index]
  se<-"Error"
  if(est.var)
  {
   if(selected=="MSCN2") se<-try(se.est(P=fit$coefficients,y=y,X=X,dist=selected,nu=fit$nu,gamma=fit$gamma),silent=TRUE)
   else se<-try(se.est(P=fit$coefficients,y=y,X=X,dist=selected,nu=fit$nu),silent=TRUE)
  } 
  if(!grepl("Error",se[1])){
  names(se)<-names(fit$coefficients)
  if(selected!="MSCN2")
  {
   if(selected!="MSN")   RVAL<- list(coefficients=fit$coefficients, se=se, nu=fit$nu, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSN", selected.model=selected)  
   else  RVAL<- list(coefficients=fit$coefficients, se=se, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSN", selected.model=selected) 
  }
  else
  {
   RVAL<- list(coefficients=fit$coefficients, se=se, nu=fit$nu, gamma=fit$gamma, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSN", selected.model=selected)  
  }
  }
  else{
  RVAL<- list(coefficients=fit$coefficients, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSN", selected.model=selected)
  RVAL$warning="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
}
 class(RVAL)<- "skewMLRM"
 RVAL$y<-y; RVAL$X<-X
 RVAL$dist<-RVAL$selected.model
 RVAL$choose.crit<-criteria
 RVAL$"function"<-"choose.MSMSN"
 RVAL
}
