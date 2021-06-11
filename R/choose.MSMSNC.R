choose.MSMSNC <-
function(y, X=NULL, max.iter=1000, prec=1e-4, est.var=TRUE, criteria="AIC",cluster=FALSE)
{
  if(cluster!=TRUE) cluster=FALSE
  if(all(criteria!=c("AIC","BIC")))
        stop("criteria must be AIC or BIC")
  if (cluster==FALSE){
  fit.msnc=try(estimate.MSNC(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mstec=try(estimate.MSTEC(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.msslec=try(estimate.MSSLEC(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  fit.mscec=try(estimate.MSCEC(y,X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)
  }
  else {
   cl0 <- parallel::detectCores()-1
   cl <- parallel::makeCluster(cl0)
  doParallel::registerDoParallel(cl);j=1
      kname=c("MSNC","MSTEC","MSSLEC","MSCEC")
  fit4 <- foreach(j = 1:4) %dopar%{
     #source('SMSNCfunctions.r')
     require(skewMLRM)
     XX=kname[j]
     estimate.dist.XX <- get(paste("estimate.", XX, sep = ""), mode = "function")
     fit4= try(estimate.dist.XX(y=y,X=X,max.iter=max.iter,prec=prec,est.var=FALSE),silent=TRUE)  
   }
   fit.msnc=fit4[[1]]
   fit.mstec=fit4[[2]]
   fit.msslec=fit4[[3]]
   fit.mscec=fit4[[4]]
  parallel::stopCluster(cl)
   }
  if(grepl("Error",fit.mstec)[1] & grepl("Error",fit.msslec)[1] & 
grepl("Error",fit.mscec)[1] & grepl("Error",fit.msnc)[1])
  {stop("estimation problem in all the models in the SMSNC class")}
  index<-c(); maxi=c()
  if(criteria=="AIC")
  {
  if(grepl("Error",fit.mstec)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mstec$AIC)}
  if(grepl("Error",fit.msslec)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.msslec$AIC)}
  if(grepl("Error",fit.mscec)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mscec$AIC)}
  if(grepl("Error",fit.msnc)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msnc$AIC)}
  }
  if(criteria=="BIC")
  {
  if(grepl("Error",fit.mstec)[1]==FALSE){ index<-c(index,1); maxi<-c(maxi, fit.mstec$BIC)}
  if(grepl("Error",fit.msslec)[1]==FALSE){ index<-c(index,2); maxi<-c(maxi, fit.msslec$BIC)}
  if(grepl("Error",fit.mscec)[1]==FALSE){ index<-c(index,3); maxi<-c(maxi, fit.mscec$BIC)}
  if(grepl("Error",fit.msnc)[1]==FALSE){ index<-c(index,4); maxi<-c(maxi, fit.msnc$BIC)}
  }
  fitted<-c("STEC","SSLEC","SCEC","SNC")[index]
  if(ncol(as.matrix(y))>1) fitted<-c("MSTEC","MSSLEC","MSCEC","MSNC")[index]
  index<-index[which.min(maxi)]
  fit<-switch(index, "1"=fit.mstec, "2"=fit.msslec, "3"=fit.mscec, "4"=fit.msnc)
  selected<- c("STEC","SSLEC","SCEC","SNC")[index]
  if(ncol(as.matrix(y))>1) selected<-c("MSTEC","MSSLEC","MSCEC","MSNC")[index]
 if(selected=="MSNC")  RVAL<- list(fitted.models=fitted, selected.model=selected, estimate=fit$estimate, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv)
 if(selected=="MSTEC" | selected=="MSSLEC")  RVAL<- list(fitted.models=fitted, selected.model=selected, estimate=fit$estimate, nu=fit$nu, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv)
 if(selected=="MSCEC")  RVAL<- list(fitted.models=fitted, selected.model=selected, estimate=fit$estimate, nu=fit$nu, gamma=fit$gamma, logLik=fit$logLik,
AIC=fit$AIC, BIC=fit$BIC, iterations=fit$iterations, conv=fit$conv)
  if(est.var)
  {
   if(selected=="MSNC") se<-try(se.est(fit$estimate[,1],y,X,dist=selected),silent=TRUE)
   if(selected=="MSTEC" | selected=="MSSLEC") se<-try(se.est(fit$estimate[,1],y,X,dist=selected,nu=fit$nu),silent=TRUE)
   if(selected=="MSCEC") se<-try(se.est(fit$estimate[,1],y,X,dist=selected,nu=fit$nu,gamma=fit$gamma),silent=TRUE)
   if(grepl("Error",se[1])) RVAL$warning="Standard errors can't be estimated: Numerical problems with the inversion of the information matrix"
   else RVAL$estimate<-cbind(RVAL$estimate,se); colnames(RVAL$estimate)<-c("estimate","s.e.")
  } 
  RVAL
}
