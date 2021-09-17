choose.models <-
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
  fit.msmn=try(choose.MSMN(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mssmn=try(choose.MSSMN(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.msmsn=try(choose.MSMSN(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.msmsnc=try(choose.MSMSNC(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  }
  else {
   cl0 <- parallel::detectCores()-1
   cl <- parallel::makeCluster(cl0)
  doParallel::registerDoParallel(cl);j=1
      kname=c("MSMN","MSSMN","MSMSN","MSMSNC")
  fit4 <- foreach(j = 1:4) %dopar%{
     #source('choosefunctions.r')
     require(skewMLRM)
     XX=kname[j]
     choose.XX <- get(paste("choose.", XX, sep = ""), mode = "function")
     fit4= try(choose.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)  
   }
   fit.msmn=fit4[[1]]
   fit.mssmn=fit4[[2]]
   fit.msmsn=fit4[[3]]
   fit.msmsnc=fit4[[4]]
  parallel::stopCluster(cl)
   }
  if(grepl("Error",fit.msmn)[1] & grepl("Error",fit.mssmn)[1] & 
     grepl("Error",fit.msmsn)[1] & grepl("Error",fit.msmsnc)[1])
  {stop("estimation problem in all the models in the MSMN class")}
  index<-c(); maxi=c()
  if(criteria=="AIC")
  {
  if(grepl("Error",fit.msmn)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.msmn$AIC)}
  if(grepl("Error",fit.mssmn)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.mssmn$AIC)}
  if(grepl("Error",fit.msmsn)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.msmsn$AIC)}
  if(grepl("Error",fit.msmsnc)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msmsnc$AIC)}
  }
  if(criteria=="BIC")
  {
  if(grepl("Error",fit.msmn)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.msmn$BIC)}
  if(grepl("Error",fit.mssmn)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.mssmn$BIC)}
  if(grepl("Error",fit.msmsn)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.msmsn$BIC)}
  if(grepl("Error",fit.msmsnc)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msmsnc$BIC)}
  }
  fitted<-c(fit.msmn$fitted.models, fit.mssmn$fitted.models,
fit.msmsn$fitted.models, fit.msmsnc$fitted.models)
  sel<-c(fit.msmn$selected.model, fit.mssmn$selected.model,
fit.msmsn$selected.model, fit.msmsnc$selected.model)[index]
  index<-index[which.min(maxi)]
  fit<-switch(index, "1"=fit.msmn, "2"=fit.mssmn, "3"=fit.msmsn, "4"=fit.msmsnc)
  selected<- sel[index]
  se<-"Error"
  if(est.var)
  {
  if(fit$dist!="MSTEC" & fit$dist!="MSSLEC" & fit$dist!="MSCEC" & fit$dist!="MSTT" & fit$dist!="MSSL2" & fit$dist!="MSCN2") se<-try(se.est(fit$coefficients,y,fit$X,dist=fit$dist),silent=TRUE)
  if(fit$dist=="MSTEC" | fit$dist=="MSSLEC" | fit$dist=="MSTT" | fit$dist=="MSSL2") se<-try(se.est(fit$coefficients,y,fit$X,dist=fit$dist,nu=fit$nu),silent=TRUE)
  if(fit$dist=="MSCEC" | fit$dist=="MSCN2") se<-try(se.est(fit$coefficients,y,fit$X,dist=fit$dist,nu=fit$nu,gamma=fit$gamma),silent=TRUE)
  } 
  if(!grepl("Error",se[1])){
  names(se)<-names(fit$coefficients)
 if (index == 4) {
  if(selected!="MSCEC")
  {
   if(selected!="MSNC")   RVAL<- list(coefficients=fit$coefficients, se=se, nu=fit$nu, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSNC", selected.model=selected)  
   if(selected=="MSNC")  RVAL<- list(coefficients=fit$coefficients, se=se, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSNC", selected.model=selected)  
  }
  if(selected=="MSCEC")   RVAL<- list(coefficients=fit$coefficients, se=se, nu=fit$nu, gamma=fit$gamma, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSNC", selected.model=selected)  
  }
  if (index == 1 | index==2) {
        RVAL<- list(coefficients=fit$coefficients, se=se, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMN", selected.model=selected)  
    }
  if (index == 3) {
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
  }
  else
  {
 if (index == 4) {
  if(selected!="MSCEC")
  {
   if(selected!="MSNC")   RVAL<- list(coefficients=fit$coefficients, nu=fit$nu, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSNC", selected.model=selected)  
   if(selected=="MSNC")  RVAL<- list(coefficients=fit$coefficients, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSNC", selected.model=selected)  
  }
  if(selected=="MSCEC")   RVAL<- list(coefficients=fit$coefficients, nu=fit$nu, gamma=fit$gamma, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSNC", selected.model=selected)  
  }
  if (index == 1 | index==2) {
        RVAL<- list(coefficients=fit$coefficients, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMN", selected.model=selected)  
    }
  if (index == 3) {
 if(selected!="MSCN2")
  {
   if(selected!="MSN")   RVAL<- list(coefficients=fit$coefficients, nu=fit$nu, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSN", selected.model=selected)  
   else  RVAL<- list(coefficients=fit$coefficients, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSN", selected.model=selected) 
  }
  else
  {
   RVAL<- list(coefficients=fit$coefficients, nu=fit$nu, gamma=fit$gamma, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, fitted.models=fitted, 
      class="MSMSN", selected.model=selected)  
  }
    }
  } 
 class(RVAL)<- "skewMLRM"
 RVAL$y<-y; RVAL$X<-X
 RVAL$dist<-RVAL$selected.model
 RVAL$choose.crit<-criteria
 RVAL$"function"<-"choose.models"
 RVAL
}
