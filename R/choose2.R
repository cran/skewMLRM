choose2 <-
function(y, X=NULL, max.iter=1000, prec=1e-4, class="MSMN", est.var=TRUE, criteria="AIC",criteria.cov="AIC", significance=0.05, cluster=FALSE)
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
 test=try(solve(MI.obs,tol=1e-100),silent=TRUE)
 se=c()
 if(is.numeric(test) & max(diag(test))<0) 
 {
 se=sqrt(-diag(test))
 }
 else  stop("Standard errors can't be estimated: Numerical problems with the inversion of the information matrix")
 se
}
 if(cluster!=TRUE) cluster=FALSE
 if (!any(class == c("MSMN","MSSMN","MSMSN","MSMSNC","ALL")))
     stop("class is not recognized")
  if (!any(criteria == c("AIC","BIC")))
     stop("criteria must be AIC or BIC")
  if (!any(criteria.cov == c("AIC","BIC","significance")))
     stop("criteria for covariates must be AIC, BIC or significance")
  if(criteria.cov=="significance")
  {
   if(!is.numeric(significance))
         stop("significance must be numeric")  
   if(significance<=0)
         stop("significance must be a positive value")  
   if(significance>0.10)
         stop("significance is usually less than or equal to 0.10")  
  }
sum.ones<-function(x) sum(x==1)
pos.ones<-function(x) which(apply(x,2,sum)==1)
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
  m=sum(mapply(sum.ones,X))/(length(X)*nrow(X[[1]]))
  q=ncol(y)
  if((ncol(X[[1]])-m)==q)
       stop("X have only the intercept term(s)")
  if (ncol(y) != nrow(X[[1]]))
        stop("y does not have the same number of columns than X")
  if (nrow(y) != length(X))
        stop("y does not have the same number of observations than X")
 if(class=="MSMN") cluster=FALSE
 if(any(class == c("MSMN","MSSMN","MSMSN","MSMSNC")))
  {
    choose.XX <- get(paste("choose.", class, sep = ""), mode = "function")
    fit=try(choose.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE, cluster=cluster),silent=TRUE)
  }
  if(class=="ALL")
  {
   fit=try(choose.models(y, X, max.iter=max.iter, prec=prec, est.var=FALSE, criteria=criteria, cluster=cluster),silent=TRUE)
  }
  if(grepl("Error",fit)[1])
  {stop(paste("estimation problem in all the models in the ",class,"class"))}
  fitted<-fit$fitted
  selected<-fit$selected
  dist<-selected
  if(criteria.cov=="significance")
  {
     object<-mbacksign(y, X=X, max.iter=max.iter, prec=prec, dist=dist, significance=significance) 
  }
  else
  {
    object=mbackcrit(y, X=X, max.iter=max.iter, prec=prec, dist=dist, criteria=criteria, est.var=FALSE, cluster=cluster)
  }
fit2<-fit
  se<-"Error";   fit<-object
  if(est.var)
  {
  if(fit$dist!="MSTEC" & fit$dist!="MSSLEC" & fit$dist!="MSCEC" & fit$dist!="MSTT" & fit$dist!="MSSL2" & fit$dist!="MSCN2") se<-try(se.est(fit$coefficients,y,fit$X,dist=fit$dist),silent=TRUE)
  if(fit$dist=="MSTEC" | fit$dist=="MSSLEC" | fit$dist=="MSTT" | fit$dist=="MSSL2") se<-try(se.est(fit$coefficients,y,fit$X,dist=fit$dist,nu=fit$nu),silent=TRUE)
  if(fit$dist=="MSCEC" | fit$dist=="MSCN2") se<-try(se.est(fit$coefficients,y,fit$X,dist=fit$dist,nu=fit$nu,gamma=fit$gamma),silent=TRUE)
  } 
  if(!grepl("Error",se[1])){
  names(se)<-names(fit$coefficients)
 if (fit$class=="MSMSNC") {
  if(fit$dist!="MSCEC")
  {
   if(fit$dist!="MSNC")   RVAL<- list(coefficients=fit$coefficients, se=se, nu=fit$nu, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
   if(fit$dist=="MSNC")  RVAL<- list(coefficients=fit$coefficients, se=se, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
  }
  if(fit$dist=="MSCEC")   RVAL<- list(coefficients=fit$coefficients, se=se, nu=fit$nu, gamma=fit$gamma, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
  }
  if (fit$class=="MSMN" | fit$class=="MSSMN") {
        RVAL<- list(coefficients=fit$coefficients, se=se, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
    }
  if (fit$class=="MSMSN") {
 if(fit$dist!="MSCN2")
  {
   if(fit$class!="MSN")   RVAL<- list(coefficients=fit$coefficients, se=se, nu=fit$nu, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
   else  RVAL<- list(coefficients=fit$coefficients, se=se, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class) 
  }
  else
  {
   RVAL<- list(coefficients=fit$coefficients, se=se, nu=fit$nu, gamma=fit$gamma, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
  }
    }
  }
  else
  {
 if (fit$class=="MSMSNC") {
  if(fit$dist!="MSCEC")
  {
   if(fit$dist!="MSNC")   RVAL<- list(coefficients=fit$coefficients, nu=fit$nu, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
   if(fit$dist=="MSNC")  RVAL<- list(coefficients=fit$coefficients, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
  }
  if(fit$dist=="MSCEC")   RVAL<- list(coefficients=fit$coefficients, nu=fit$nu, gamma=fit$gamma, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
  }
  if (fit$class=="MSMN" | fit$class=="MSSMN") {
        RVAL<- list(coefficients=fit$coefficients, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
    }
  if (fit$class=="MSMSN") {
 if(fit$dist!="MSCN2")
  {
   if(fit$dist!="MSN")   RVAL<- list(coefficients=fit$coefficients, nu=fit$nu, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
   else  RVAL<- list(coefficients=fit$coefficients, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class) 
  }
  else
  {
   RVAL<- list(coefficients=fit$coefficients, nu=fit$nu, gamma=fit$gamma, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv, dist=fit$dist, class=fit$class)  
  }
    }
  } 
 class(RVAL)<- "skewMLRM"
if(!is.null(object$"eliminated"))
   {
     RVAL$"eliminated"=object$"eliminated"
   }
   if(!is.null(object$"warnings"))
   {
     RVAL$warnings=object$"warnings"
   }
 RVAL$"function"<-"choose2"
 RVAL$choose.crit<-criteria
 RVAL$choose.crit.cov<-criteria.cov
 if(RVAL$choose.crit.cov=="significance") RVAL$sign<-significance
 RVAL$y<-y; RVAL$X<-X
 RVAL$"function"<-"choose2"
 RVAL$dist<-object$dist
 RVAL$class<-object$class
 RVAL$comment <- fit2$comment
 if(!is.null(fit$eliminated)) RVAL$eliminated <- fit2$eliminated
 RVAL$fitted.models<-fit2$fitted.models
if(fit2$class=="MSMSN" | fit2$class=="MSMSNC")
{
  if(fit2$dist=="MSTT" | fit2$dist=="MSSL2" | fit2$dist=="MSTEC" | fit2$dist=="MSSLEC")
  {RVAL$nu=fit2$nu}
  if(fit2$dist=="MSCN2" | fit2$dist=="MSCEC")
  {RVAL$nu=fit2$nu; RVAL$gamma=fit2$gamma}
}
RVAL$selected.model=fit2$selected.model
RVAL$fitted.classes=class
 RVAL
}
